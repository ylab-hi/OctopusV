from pathlib import Path

import typer

from octopusv.converter.base import get_alt_chrom_pos
from octopusv.converter.bnd2dup import BND_to_DUP_Converter
from octopusv.converter.bnd2inv import BND_to_INV_Converter
from octopusv.converter.bnd2tra_forward import BND_to_TRA_Forward_Converter
from octopusv.converter.bnd2tra_reverse import BND_to_TRA_Reverse_Converter
from octopusv.converter.mpi2tra import MatePairIndependentToTRAConverter
from octopusv.converter.mpm2tra import MatePairMergeToTRAConverter
from octopusv.converter.mprtra2tra import MatePairReciprocalTranslocationToTRAConverter
from octopusv.converter.nobnd import NonBNDConverter
from octopusv.converter.snmd_dndpi2tra import SpecialNoMateDiffBNDPairIndependentToTRAConverter
from octopusv.converter.snmd_dndpr_tra2tra import SpecialNoMateDiffBNDPairReciprocalTranslocationToTRAConverter
from octopusv.converter.stra2tra import SingleTRAToTRAConverter
from octopusv.transformer.mp_bnd import MatePairBNDTransformer
from octopusv.transformer.no_bnd import NonBNDTransformer
from octopusv.transformer.same_chr_dnd import SameChrBNDTransformer
from octopusv.transformer.snmd_bndp import SpecialNoMateDiffBNDPairTransformer
from octopusv.transformer.stra import SingleTRATransformer
from octopusv.utils.parser import parse_vcf
from octopusv.utils.svcf_utils import write_sv_vcf


def correct(
    input_vcf: Path = typer.Argument(
        ...,
        exists=True,
        dir_okay=False,
        resolve_path=True,
    ),
    output: Path = typer.Argument(..., dir_okay=False, resolve_path=True),
    pos_tolerance: int = typer.Option(
        3,
        "--pos-tolerance",
        "-pt",
        help="Position tolerance for identifying mate BND events, default=3, recommend not to set larger than 5",
    ),
):
    """Correct SV events."""
    # Parse the input VCF file
    # non_bnd_events means DEL, INV, INS, DUP
    contig_lines, same_chr_bnd_events, diff_chr_bnd_events, non_bnd_events = parse_vcf(
        input_vcf,
    )

    # Extract mate BND and no mate events, they are all with different chromosomes
    mate_bnd_pairs = find_mate_bnd_events(
        diff_chr_bnd_events,
        pos_tolerance=pos_tolerance,
    )

    no_mate_events = find_no_mate_events(
        diff_chr_bnd_events,
        pos_tolerance=pos_tolerance,
    )

    # Further classify no_mate_events into special_no_mate_diff_bnd_pair and other_single_TRA
    (
        special_no_mate_diff_bnd_pairs,
        other_single_TRA,
    ) = find_special_no_mate_diff_bnd_pair_and_other_single_tra(
        no_mate_events,
        pos_tolerance=pos_tolerance,
    )

    # Initialize the EventTransformer with a list of transform strategies for each type of events
    # One transformer and multiple-strategies
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
    non_bnd_transformer = NonBNDTransformer(
        [NonBNDConverter()],
    )  # Assuming non-BND events are not to be transformed for now

    # Apply all transformation strategies to the events
    same_chr_bnd_transformed_events = same_chr_bnd_transformer.apply_transforms(
        same_chr_bnd_events,
    )
    mate_pair_transformed_events = mate_bnd_pair_transformer.apply_transforms(
        mate_bnd_pairs,
    )
    special_no_mate_diff_bnd_pair_transformed_events = special_no_mate_diff_bnd_pair_transformer.apply_transforms(
        special_no_mate_diff_bnd_pairs,
    )
    single_TRA_transformed_events = single_TRA_transformer.apply_transforms(
        other_single_TRA,
    )
    non_bnd_transformed_events = non_bnd_transformer.apply_transforms(non_bnd_events)

    # Merge all transformed events
    all_transformed_events = (
        same_chr_bnd_transformed_events
        + mate_pair_transformed_events
        + special_no_mate_diff_bnd_pair_transformed_events
        + single_TRA_transformed_events
        + non_bnd_transformed_events
    )

    # Write the transformed events to the output file
    write_sv_vcf(contig_lines, all_transformed_events, output)

# ==============================
# The following are auxiliary Functions.
# ==============================
def find_mate_bnd_events(events, pos_tolerance=3):
    """Find mate BND events."""
    event_dict = {}
    mate_bnd_pairs = []

    for event in events:
        chrom_alt, pos_alt = get_alt_chrom_pos(event.alt)
        key_for_searching = (event.chrom, event.pos, chrom_alt, pos_alt)

        # Generate all possible reverse keys
        possible_reverse_keys = [
            (chrom_alt, pos_alt + i, event.chrom, event.pos + j)
            for i in range(-pos_tolerance, pos_tolerance + 1)
            for j in range(-pos_tolerance, pos_tolerance + 1)
        ]

        for reverse_key in possible_reverse_keys:
            if reverse_key in event_dict:
                mate_bnd_pairs.append((event_dict.pop(reverse_key), event))
                break
        else:
            event_dict[key_for_searching] = event
            # It will override for no_mate_events due to the same key, but does matter.
            # We only output mate_bnd_pairs here

    return mate_bnd_pairs


def find_no_mate_events(events, pos_tolerance=3):
    """Find no mate events. The main idea is to store events in list to avoid overriding in event_dict."""
    event_dict = {}
    no_mate_events = []
    mate_bnd_pairs = []

    for event in events:
        chrom_alt, pos_alt = get_alt_chrom_pos(event.alt)
        key = (event.chrom, event.pos, chrom_alt, pos_alt)

        # Generate all possible reverse keys
        possible_reverse_keys = [
            (chrom_alt, pos_alt + i, event.chrom, event.pos + j)
            for i in range(-pos_tolerance, pos_tolerance + 1)
            for j in range(-pos_tolerance, pos_tolerance + 1)
        ]

        mate_found = False  # A flag to indicate if a matching mate has been found for the event.

        # Check each reverse key to see if a mate exists in event_dict.
        for reverse_key in possible_reverse_keys:
            if reverse_key in event_dict:   # Means got a mate.
                # Now iterating through all events associated with the reverse_key
                for mate_event in event_dict[reverse_key]:
                    # If a mate is found, add the pair to the mate_bnd_pairs list.
                    mate_bnd_pairs.append((mate_event, event))
                    # Remove the found mate event from the dictionary.
                    event_dict[reverse_key].remove(mate_event)
                    mate_found = True
                    break  # Stop searching once a mate is found.
                # If a mate was found and no events are left under the reverse key, remove the key.
                if mate_found:
                    if len(event_dict[reverse_key]) == 0:
                        # Remove the key if no events are left
                        del event_dict[reverse_key]
                    break

        if not mate_found:
            if key not in event_dict:  # Create a new list if key doesn't exist.
                event_dict[key] = []
            event_dict[key].append(event)  # # Add the event to the list under its key.

    # Flattening the list of events from all keys
    no_mate_events = [single_event for event_list in event_dict.values() for single_event in event_list]

    return no_mate_events


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
            if reverse_key in event_dict:  # If key in dict.
                special_no_mate_diff_bnd_pair.append(
                    (event_dict.pop(reverse_key), event),
                )  # event_dict.pop(reverse_key) will delete mate events from event_dic and output
                mate_found = True
                break

        if not mate_found:
            event_dict[key] = event

    other_single_TRA = list(event_dict.values())

    return special_no_mate_diff_bnd_pair, other_single_TRA
