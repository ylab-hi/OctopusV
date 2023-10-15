import argparse


def main():
    parser = argparse.ArgumentParser(description="SVmerger tool")
    parser.add_argument("mode", type=str, help="Mode of operation")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input VCF file")
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Output converted VCF file",
    )
    parser.add_argument(
        "-p",
        "--pos-tolerance",
        type=int,
        default=3,
        help="Position tolerance for identifying mate BND events, default=3, recommend not to set larger than 5",
    )

    args = parser.parse_args()

    if args.mode == "convert":
        # Parse the input VCF file
        headers, same_chr_bnd_events, diff_chr_bnd_events, non_bnd_events = parse_vcf(
            args.input,
        )  # non_bnd_events means DEL, INV, INS, DUP

        # Extract mate BND and no mate events, they are all with different chromosomes
        mate_bnd_pairs, no_mate_events = find_mate_bnd_and_no_mate_events(
            diff_chr_bnd_events,
            pos_tolerance=args.pos_tolerance,
        )

        # Further classify no_mate_events into special_no_mate_diff_bnd_pair and other_single_TRA
        (
            special_no_mate_diff_bnd_pairs,
            other_single_TRA,
        ) = find_special_no_mate_diff_bnd_pair_and_other_single_tra(
            no_mate_events,
            pos_tolerance=args.pos_tolerance,
        )

        # Initialize the EventTransformer with a list of transform strategies for each type of events: one transformer and multiple-strategies
        same_chr_bnd_transformer = SameChrBNDTransformer(
            [
                BND_to_INV_Converter(),
                BND_to_DUP_Converter(),
                BND_to_TRA_Forward_Converter(),
                BND_to_TRA_Reverse_Converter(),
            ],
        )
        mate_bnd_pair_transformer = MatePairBNDTransformer(
            [
                MatePairReciprocalTranslocationToTRAConverter(),
                MatePairIndependentToTRAConverter(),
                MatePairMergeToTRAConverter(),
            ],
        )
        special_no_mate_diff_bnd_pair_transformer = SpecialNoMateDiffBNDPairTransformer(
            [
                SpecialNoMateDiffBNDPairReciprocalTranslocationToTRAConverter(),
                SpecialNoMateDiffBNDPairIndependentToTRAConverter(),
            ],
        )
        single_TRA_transformer = SingleTRATransformer([SingleTRAToTRAConverter()])
        non_bnd_transformer = EventTransformer(
            [],
        )  # Assuming non-BND events are not to be transformed

        # Apply all transformation strategies to the events
        same_chr_bnd_transformed_events = same_chr_bnd_transformer.apply_transforms(
            same_chr_bnd_events,
        )
        mate_pair_transformed_events = mate_bnd_pair_transformer.apply_transforms(
            mate_bnd_pairs,
        )
        special_no_mate_diff_bnd_pair_transformed_events = (
            special_no_mate_diff_bnd_pair_transformer.apply_transforms(
                special_no_mate_diff_bnd_pairs,
            )
        )
        single_TRA_transformed_events = single_TRA_transformer.apply_transforms(
            other_single_TRA,
        )
        non_bnd_transformer.apply_transforms(non_bnd_events)

        # Merge all transformed events
        all_transformed_events = (
            same_chr_bnd_transformed_events
            + mate_pair_transformed_events
            + special_no_mate_diff_bnd_pair_transformed_events
            + single_TRA_transformed_events
            + non_bnd_events
        )

        # Write the transformed events to the output file
        same_chr_bnd_transformer.write_vcf(headers, all_transformed_events, args.output)
    else:
        print(f"Unknown mode: {args.mode}")
        parser.print_help()


# ==============================
# The following are auxiliary Functions.
# ==============================


def find_mate_bnd_and_no_mate_events(events, pos_tolerance=3):
    """Extract mate BND and no mate events from diff_chr_bnd_events."""
    event_dict = {}
    mate_bnd_pairs = []
    no_mate_events = []

    for event in events:
        chrom_alt, pos_alt = get_alt_chrom_pos(event.alt)
        key = (event.chrom, event.pos, chrom_alt, pos_alt)

        # Generate all possible reverse keys
        possible_reverse_keys = [
            (chrom_alt, pos_alt + i, event.chrom, event.pos + j)
            for i in range(-pos_tolerance, pos_tolerance + 1)
            for j in range(-pos_tolerance, pos_tolerance + 1)
        ]

        mate_found = False  # This is a flag
        for reverse_key in possible_reverse_keys:
            if reverse_key in event_dict:
                mate_bnd_pairs.append(
                    (event_dict.pop(reverse_key), event),
                )  # event_dict.pop(reverse_key) will delete mate events from event_dic and output
                mate_found = True
                break

        if not mate_found:
            event_dict[key] = event

    no_mate_events = list(event_dict.values())

    return mate_bnd_pairs, no_mate_events


def find_special_no_mate_diff_bnd_pair_and_other_single_tra(
    events,
    pos_tolerance=3,
):  # Process no_mate_events
    """Extract special no mate diff bnd pair and other single TRA events from no_mate_events."""
    event_dict = {}
    special_no_mate_diff_bnd_pair = []
    other_single_TRA = []

    for event in events:
        chrom_alt, pos_alt = get_alt_chrom_pos(event.alt)
        key = (event.chrom, event.pos, chrom_alt, pos_alt)

        # Generate all possible reverse keys
        possible_reverse_keys = [
            (event.chrom, event.pos + i, chrom_alt, pos_alt + j)
            for i in range(-pos_tolerance, pos_tolerance + 1)
            for j in range(-pos_tolerance, pos_tolerance + 1)
        ]

        mate_found = False  # This is a flag
        for reverse_key in possible_reverse_keys:
            if reverse_key in event_dict:
                special_no_mate_diff_bnd_pair.append(
                    (event_dict.pop(reverse_key), event),
                )  # event_dict.pop(reverse_key) will delete mate events from event_dic and output
                mate_found = True
                break

        if not mate_found:
            event_dict[key] = event

    other_single_TRA = list(event_dict.values())

    return special_no_mate_diff_bnd_pair, other_single_TRA


if __name__ == "__main__":
    main()

# Each software needs to calibrate the calculation method for SVLEN internally.
# After each software's conversion, redundancy must be removed internally.
# My standard format: Use <INV> for INV, use <DUP> for DUP, and for INS, DEL, it's best to keep the actual sequence.
# Now, only the non_bnd_events haven't been processed. They are actually more standard events like INS, INV, DEL, and DUP.
