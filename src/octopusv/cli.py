# -*- coding: utf-8 -*-

import argparse
import re
from natsort import natsorted


def main():
    parser = argparse.ArgumentParser(description="SVmerger tool")
    parser.add_argument("mode", type=str, help="Mode of operation")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input VCF file")
    parser.add_argument(
        "-o", "--output", type=str, required=True, help="Output converted VCF file"
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
            args.input
        )  # non_bnd_events means DEL, INV, INS, DUP

        # Extract mate BND and no mate events, they are all with different chromosomes
        mate_bnd_pairs, no_mate_events = find_mate_bnd_and_no_mate_events(
            diff_chr_bnd_events, pos_tolerance=args.pos_tolerance
        )

        # Further classify no_mate_events into special_no_mate_diff_bnd_pair and other_single_TRA
        (
            special_no_mate_diff_bnd_pairs,
            other_single_TRA,
        ) = find_special_no_mate_diff_bnd_pair_and_other_single_TRA(
            no_mate_events, pos_tolerance=args.pos_tolerance
        )

        # Initialize the EventTransformer with a list of transform strategies for each type of events: one transformer and multiple-strategies
        same_chr_bnd_transformer = SameChrBNDTransformer(
            [
                BND_to_INV_Converter(),
                BND_to_DUP_Converter(),
                BND_to_TRA_Forward_Converter(),
                BND_to_TRA_Reverse_Converter(),
            ]
        )
        mate_bnd_pair_transformer = MatePairBNDTransformer(
            [
                MatePairReciprocalTranslocationToTRAConverter(),
                MatePairIndependentToTRAConverter(),
                MatePairMergeToTRAConverter(),
            ]
        )
        special_no_mate_diff_bnd_pair_transformer = SpecialNoMateDiffBNDPairTransformer(
            [
                SpecialNoMateDiffBNDPairReciprocalTranslocationToTRAConverter(),
                SpecialNoMateDiffBNDPairIndependentToTRAConverter(),
            ]
        )
        single_TRA_transformer = SingleTRATransformer([SingleTRAToTRAConverter()])
        non_bnd_transformer = EventTransformer(
            []
        )  # Assuming non-BND events are not to be transformed

        # Apply all transformation strategies to the events
        same_chr_bnd_transformed_events = same_chr_bnd_transformer.apply_transforms(
            same_chr_bnd_events
        )
        mate_pair_transformed_events = mate_bnd_pair_transformer.apply_transforms(
            mate_bnd_pairs
        )
        special_no_mate_diff_bnd_pair_transformed_events = (
            special_no_mate_diff_bnd_pair_transformer.apply_transforms(
                special_no_mate_diff_bnd_pairs
            )
        )
        single_TRA_transformed_events = single_TRA_transformer.apply_transforms(
            other_single_TRA
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


def check_vcf_format(vcf_file_path):
    """
    Check the format of a VCF file.
    Raise an error and exit if the format is incorrect.
    """
    with open(vcf_file_path) as f:
        lines = f.readlines()

    # Check if there is at least one header line
    if not any(line.startswith("#") for line in lines):
        print(
            "ERROR: Invalid VCF format. The file should contain at least one header line starting with '#'."
        )
        exit(1)

    for line in lines:
        if line.startswith("#"):
            continue  # Skip header lines

        # Check for space in lines
        if " " in line:
            print(
                "ERROR: Invalid VCF format. Non-header lines should not contain spaces."
            )
            exit(1)

        fields = line.split("\t")

        # Check the number of columns
        if len(fields) < 10:
            print(
                f"ERROR: Invalid VCF format. Expected at least 10 fields, but got {len(fields)}"
            )
            exit(1)

        # Check that the position is a number
        try:
            int(fields[1])
        except ValueError:
            print(
                f"ERROR: Invalid VCF format. Position (field 2) should be a number, but got {fields[1]}"
            )
            exit(1)

        # Check that the quality score is a number or '.'
        if fields[5] != ".":
            try:
                float(fields[5])
            except ValueError:
                print(
                    f"ERROR: Invalid VCF format. Quality score (field 6) should be a number or '.', but got {fields[5]}"
                )
                exit(1)


def is_same_chr_bnd(event):
    """
    Check if the POS and ALT of an event are on the same chromosome.
    """
    if event.is_BND():
        split_result = re.split(r"[\[\]:]", event.alt)
        if len(split_result) != 4:
            print(
                f"Unexpected ALT format, it should be something like N]chr10:69650962]: {split_result}"
            )
        else:
            chrom_alt, _ = split_result[1:3]
            return event.chrom == chrom_alt

    return False  # For non-BND, we won't categorize them as same_chr_bnd or diff_chr_bnd events


def get_BND_pattern(alt):
    # Get the pattern of BND: t[p[, t]p], ]p]t, [p[t
    if alt[0] in "ATCGN" and alt[1] == "[":
        return "t[p["
    elif alt[0] in "ATCGN" and alt[1] == "]":
        return "t]p]"
    elif alt[-1] in "ATCGN" and alt[-2] == "[":
        return "[p[t"
    elif alt[-1] in "ATCGN" and alt[-2] == "]":
        return "]p]t"
    else:
        return None


def get_alt_chrom_pos(alt):
    """
    Extract chromosome and position from alt field in VCF file.
    """
    split_result = re.split(r"[\[\]:]", alt)
    if len(split_result) != 4:
        print(
            f"Unexpected ALT format, it should be something like N]chr10:69650962]: {split_result}"
        )
        return None, None
    else:
        chrom_alt, pos_alt = split_result[1:3]
        return chrom_alt, int(pos_alt)  # Convert pos_alt to integer.


def find_mate_bnd_and_no_mate_events(events, pos_tolerance=3):
    """
    Extract mate BND and no mate events from diff_chr_bnd_events
    """
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
                    (event_dict.pop(reverse_key), event)
                )  # event_dict.pop(reverse_key) will delete mate events from event_dic and output
                mate_found = True
                break

        if not mate_found:
            event_dict[key] = event

    no_mate_events = list(event_dict.values())

    return mate_bnd_pairs, no_mate_events


def is_MatePairReciprocalTranslocation(event1, event2):
    # Define the set of qualified pairings of BND patterns, used by MatePairReciprocalTranslocationToTRAConverter
    qualified_pairings = [
        {"]p]t", "]p]t"},
        {"t[p[", "t[p["},
        {"t]p]", "[p[t"},
    ]
    # Extract the patterns from each event
    event1_pattern = get_BND_pattern(event1.alt)
    event2_pattern = get_BND_pattern(event2.alt)

    # Check if the patterns match one of the qualified pairings
    for pairing in qualified_pairings:
        if {event1_pattern, event2_pattern} == pairing:
            return True

    return False


def is_same_bnd_event(event1, event2):
    # Define whether the BND represent the same TRA events, used by MatePairMergeToTRAConverter
    qualified_pairings = [{"]p]t", "t[p["}]
    event1_pattern = get_BND_pattern(event1.alt)
    event2_pattern = get_BND_pattern(event2.alt)

    for pairing in qualified_pairings:
        if {event1_pattern, event2_pattern} == pairing:
            return True

    return False


def compare_chromosomes(event1, event2):  # Used by MatePairMergeToTRAConverter
    """
    Determine if the chromosome of event1 is smaller than the chromosome of event2 using natural sort order.

    Returns:
    A boolean value. True if the chromosome of event1 is smaller or equal to the chromosome of event2, False otherwise.
    """
    original_order = [event1.chrom, event2.chrom]

    # Get the order when the chromosomes are naturally sorted
    sorted_order = natsorted(original_order)

    # it means the chromosome of event1 is smaller or equal to the chromosome of event2
    return original_order == sorted_order


def is_independent_bnd_event(event1, event2):  # used by MatePairMergeToTRAConverter
    """
    Determine if two mate_pair_bnd events are independent events.
    """
    qualified_pairings = [
        {"t]p]", "]p]t"},
        {"t]p]", "t]p]"},
        {"[p[t", "[p[t"},
        {"t[p[", "t]p]"},
        {"t[p[", "[p[t"},
        {"]p]t", "[p[t"},
    ]
    event1_pattern = get_BND_pattern(event1.alt)
    event2_pattern = get_BND_pattern(event2.alt)

    for pairing in qualified_pairings:
        if {event1_pattern, event2_pattern} == pairing:
            return True

    return False


def find_special_no_mate_diff_bnd_pair_and_other_single_TRA(
    events, pos_tolerance=3
):  # Process no_mate_events
    """
    Extract special no mate diff bnd pair and other single TRA events from no_mate_events.
    """
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
                    (event_dict.pop(reverse_key), event)
                )  # event_dict.pop(reverse_key) will delete mate events from event_dic and output
                mate_found = True
                break

        if not mate_found:
            event_dict[key] = event

    other_single_TRA = list(event_dict.values())

    return special_no_mate_diff_bnd_pair, other_single_TRA


def is_SpecialNoMateDiffBndPairReciprocalTranslocation(
    event1, event2
):  # Used by SpecialNoMateDiffBNDPairReciprocalTranslocationToTRAConverter
    qualified_pairings = [{"t[p[", "]p]t"}, {"t]p]", "[p[t"}]
    # Extract the patterns from each event
    event1_pattern = get_BND_pattern(event1.alt)
    event2_pattern = get_BND_pattern(event2.alt)

    # Check if the patterns match one of the qualified pairings
    for pairing in qualified_pairings:
        if {event1_pattern, event2_pattern} == pairing:
            return True

    return False


def is_independent_special_bnd_event(
    event1, event2
):  # Used by SpecialNoMateDiffBNDPairIndependentToTRAConverter
    """
    Determine if two special_no_mate_diff_bnd_pair events are independent events.
    """
    qualified_pairings = [
        {"t[p[", "t]p]"},
        {"t[p[", "[p[t"},
        {"t]p]", "]p]t"},
        {"]p]t", "[p[t"},
    ]
    event1_pattern = get_BND_pattern(event1.alt)
    event2_pattern = get_BND_pattern(event2.alt)

    for pairing in qualified_pairings:
        if {event1_pattern, event2_pattern} == pairing:
            return True

    return False


def parse_vcf(vcf_file_path):
    """
    Parse VCF file into lists of SVEvent objects based on their type (same chromosome BND, different chromosome BND, non-BND) and return headers.
    """
    check_vcf_format(vcf_file_path)  # Check the format first
    same_chr_bnd_events = []
    diff_chr_bnd_events = []
    non_bnd_events = []
    headers = []

    with open(vcf_file_path) as f:
        for line in f:
            if line.startswith("#"):
                headers.append(line.strip())
                continue  # Skip header lines

            fields = line.strip().split("\t")
            event = SVEvent(*fields)  # Unpack fields and send to SVEvent class

            if event.is_BND():
                if is_same_chr_bnd(event):  # Check if the event is same chromosome
                    same_chr_bnd_events.append(event)
                else:
                    diff_chr_bnd_events.append(event)  # Different chromosome
            else:
                non_bnd_events.append(event)  # Non-BND events

    return headers, same_chr_bnd_events, diff_chr_bnd_events, non_bnd_events


class SVEvent:  # A class to represent each SV event
    def __init__(self, chrom, pos, id, ref, alt, qual, filter, info, format, sample):
        self.chrom = chrom
        self.pos = int(pos)
        self.id = id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter
        self.info = self._parse_info(info)  # self.info will be a dict
        self.format = format
        self.sample = sample

    def _parse_info(self, info):  # Parse the info field according to your requirement
        info_dict = {}
        for item in info.split(";"):
            parts = item.split("=")
            if len(parts) == 2:  # Just in case that the INFO doesn't contain = symbol
                key, value = parts
                info_dict[key] = value
            else:
                info_dict[parts[0]] = None
        return info_dict

    def is_duplication(self):
        # Whether is duplication
        return self.info["SVTYPE"] == "DUP"

    def is_inversion(self):
        # Whether is inversion
        return self.info["SVTYPE"] == "INV"

    def is_insertion(self):
        #  Whether is insertion
        return self.info["SVTYPE"] == "INS"

    def is_TRA(self):
        #  Whether is TRA
        return self.info["SVTYPE"] == "TRA"

    def is_BND(self):
        #  Whether is BND
        return self.info["SVTYPE"] == "BND"

    def __str__(self):
        # enable you to output object as vcf string format
        # Define your desired order
        ordered_keys = ["SVTYPE", "END", "SVLEN", "SUPPORT", "COVERAGE", "STRAND"]
        info_items = []

        # First add the ordered keys
        for key in ordered_keys:
            if key in self.info:
                info_items.append(f"{key}={self.info[key]}")

        # Then add the remaining keys
        for key in self.info:
            if key not in ordered_keys:
                info_items.append(f"{key}={self.info[key]}")

        info_str = ";".join(info_items)  # SVTYPE=INV;END=69650961;SVLEN=1835141

        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            self.chrom,
            self.pos,
            self.id,
            self.ref,
            self.alt,
            self.qual,
            self.filter,
            info_str,
            self.format,
            self.sample,
        )


class Converter:
    """
    This is an abstract base class for all converter classes. It provides a common interface for all converters.
    The `convert` method is a placeholder that needs to be overridden in each concrete converter class.
    """

    def convert(self, event):
        """
        This method should be overridden in a subclass.
        """
        raise NotImplementedError


# ==============================
# The following four strategies converts a single event and does not need to return it,
# because the changes are made directly to the mutable event object.
# They are used by SameChrBNDTransformer
# ==============================


class BND_to_INV_Converter(Converter):
    """
    This class inherits from the `Converter` base class and implements the conversion logic for BND to INV conversion.
    """

    def convert(self, event):  # Each converter has a convert function
        try:
            if event.is_BND():
                pattern = get_BND_pattern(event.alt)
                if pattern == "t]p]" or pattern == "[p[t":
                    chrom_alt, pos_alt = get_alt_chrom_pos(event.alt)
                    if chrom_alt is None:
                        print("Failed to get ALT chrom and pos")
                    else:
                        if (
                            event.chrom == chrom_alt
                        ):  # Do this only when same chromosome
                            if pattern == "t]p]" and event.pos < pos_alt:
                                end = pos_alt
                                svlen = abs(event.pos - pos_alt)
                                self.convert_to_INV(
                                    event, end, svlen
                                )  # self will be passed to convert_to_INV method.
                            elif pattern == "[p[t" and event.pos > pos_alt:
                                end = event.pos
                                event.pos = pos_alt
                                svlen = abs(event.pos - end)
                                self.convert_to_INV(event, end, svlen)

        except Exception as e:
            print("Failed to convert BND to Inversion: ", e)

    def convert_to_INV(self, event, end, svlen):
        # Detail logic changing BND event to INV event
        event.alt = "<INV>"
        event.info["SVTYPE"] = "INV"
        event.info["END"] = end
        event.info["SVLEN"] = svlen


class BND_to_DUP_Converter(Converter):
    """
    This class inherits from the `Converter` base class and implements the conversion logic for BND to DUP conversion.
    """

    def convert(self, event):
        try:
            if event.is_BND():
                pattern = get_BND_pattern(event.alt)
                if pattern in ["]p]t", "t[p["]:
                    chrom_alt, pos_alt = get_alt_chrom_pos(event.alt)
                    if chrom_alt is None:
                        print("Failed to get ALT chrom and pos")
                    else:
                        if event.chrom == chrom_alt:
                            if pattern == "]p]t" and event.pos < pos_alt:
                                end = pos_alt
                                svlen = abs(end - event.pos)
                                self.convert_to_DUP(event, end, svlen)
                            elif pattern == "t[p[" and event.pos > pos_alt:
                                end = event.pos
                                svlen = abs(event.pos - pos_alt)
                                event.pos = pos_alt
                                self.convert_to_DUP(event, end, svlen)
        except Exception as e:
            print("Failed to convert BND to Duplication: ", e)

    def convert_to_DUP(self, event, end, svlen):
        event.alt = "<DUP>"
        event.info["SVTYPE"] = "DUP"
        event.info["END"] = end
        event.info["SVLEN"] = svlen


class BND_to_TRA_Forward_Converter(Converter):
    """
    This class inherits from the `Converter` base class and implements
    the conversion logic for BND to TRA Forward translocation cut-paste.
    """

    def convert(self, event):
        try:
            if event.is_BND():
                pattern = get_BND_pattern(event.alt)
                if pattern in ["t[p[", "]p]t"]:
                    chrom_alt, pos_alt = get_alt_chrom_pos(event.alt)
                    if chrom_alt is None:
                        print("Failed to get ALT chrom and pos")
                    else:
                        if event.chrom == chrom_alt:
                            if pattern == "t[p[" and event.pos < pos_alt:
                                end = pos_alt
                                self.convert_to_TRA_forward(event, end)
                            elif pattern == "]p]t" and event.pos > pos_alt:
                                end = pos_alt
                                self.convert_to_TRA_forward(event, end)
        except Exception as e:
            print("Failed to convert BND to Translocation: ", e)

    def convert_to_TRA_forward(self, event, end):
        event.info["SVTYPE"] = "TRA"
        event.info["END"] = end
        event.info["SVLEN"] = 0
        event.info["CHR2"] = event.chrom


class BND_to_TRA_Reverse_Converter(Converter):
    """
    This class inherits from the `Converter` base class and implements the conversion logic for BND to TRA reverse conversion.
    """

    def convert(self, event):
        try:
            if event.is_BND():
                pattern = get_BND_pattern(event.alt)
                if pattern in ["[p[t", "t]p]"]:
                    chrom_alt, pos_alt = get_alt_chrom_pos(event.alt)
                    if chrom_alt is None:
                        print("Failed to get ALT chrom and pos")
                    else:
                        if event.chrom == chrom_alt:
                            if pattern == "t]p]" and event.pos > pos_alt:
                                end = pos_alt
                                self.convert_to_TRA_reverse(event, end)
                            elif pattern == "[p[t" and event.pos < pos_alt:
                                end = pos_alt
                                self.convert_to_TRA_reverse(event, end)
        except Exception as e:
            print("Failed to convert BND to Translocation: ", e)

    def convert_to_TRA_reverse(self, event, end):
        event.info["SVTYPE"] = "TRA"
        event.info["END"] = end
        event.info["SVLEN"] = 0
        event.info["CHR2"] = event.chrom


# You can add more converter classes here...
# ==============================
# The following strategies convert Mate pair bnd events to reciprocal translocation, independent translocation,
# or merge the same events.
# ==============================


class MatePairReciprocalTranslocationToTRAConverter(Converter):
    """
    This class inherits from the `Converter` base class and implements
    the conversion logic for reciprocal translocation.
    """

    def convert(self, pair):
        event1, event2 = pair
        # Check if this pair satisfies the criteria for reciprocal translocation
        if is_MatePairReciprocalTranslocation(event1, event2):
            # Convert events
            self.convert_to_TRA(event1, event2)
            return [event1, event2]  # Return a list of transformed events
        else:
            return []  # If the pair doesn't satisfy the criteria, return an empty list

    def convert_to_TRA(self, event1, event2):
        # Convert a pair of events to reciprocal translocation
        # Modify event1
        event1.info["SVTYPE"] = "TRA"
        chrom_alt, pos_alt = get_alt_chrom_pos(event1.alt)
        event1.info["END"] = pos_alt
        event1.info["CHR2"] = chrom_alt
        event1.info["RTID"] = event2.id
        event1.info["SVLEN"] = 0
        # Modify event2
        event2.info["SVTYPE"] = "TRA"
        chrom_alt, pos_alt = get_alt_chrom_pos(event2.alt)
        event2.info["END"] = pos_alt
        event2.info["CHR2"] = chrom_alt
        event2.info["RTID"] = event1.id
        event2.info["SVLEN"] = 0
        # No need to return anything as the states of event1 and event2 are modified in-place


class MatePairMergeToTRAConverter(Converter):
    """Class to transform a pair of identical BND events into a single TRA event."""

    def convert(self, pair):
        event1, event2 = pair
        # Check if the pair of events has the same BND pattern
        if is_same_bnd_event(event1, event2):
            # Determine which event should be retained
            if compare_chromosomes(event1, event2):  # If the first event is "smaller"
                retained_event = event1
            else:
                retained_event = event2

            # Modify the SVTYPE to TRA
            retained_event.info["SVTYPE"] = "TRA"
            retained_event.info["SVLEN"] = 0
            chrom_alt, pos_alt = get_alt_chrom_pos(retained_event.alt)
            retained_event.info["END"] = pos_alt
            retained_event.info["CHR2"] = chrom_alt

            # Return the retained event only
            return [retained_event]
        else:
            # If the pair of events does not have the same BND pattern, return an empty list
            return []


class MatePairIndependentToTRAConverter(Converter):
    """
    the conversion logic for independent translocation.
    """

    def convert(self, pair):
        event1, event2 = pair
        # Check if this pair satisfies the criteria for independent translocation
        if is_independent_bnd_event(event1, event2):
            # Convert events
            self.convert_to_TRA(event1, event2)
            return [event1, event2]  # Return a list of transformed events
        else:
            return []  # If the pair doesn't satisfy the criteria, return an empty list

    def convert_to_TRA(self, event1, event2):
        # Convert a pair of independent events to TRA
        # Modify event1
        event1.info["SVTYPE"] = "TRA"
        event1.info["SVLEN"] = 0
        chrom_alt, pos_alt = get_alt_chrom_pos(event1.alt)
        event1.info["END"] = pos_alt
        event1.info["CHR2"] = chrom_alt
        # Modify event2
        event2.info["SVTYPE"] = "TRA"
        event2.info["SVLEN"] = 0
        chrom_alt, pos_alt = get_alt_chrom_pos(event2.alt)
        event2.info["END"] = pos_alt
        event2.info["CHR2"] = chrom_alt
        # No need to return anything as the states of event1 and event2 are modified in-place


# ==============================
# The following strategies convert special_no_mate_diff_bnd_pair to reciprocal translocation, independent translocation
# ==============================
class SpecialNoMateDiffBNDPairReciprocalTranslocationToTRAConverter(Converter):
    """
    This class inherits from the `Converter` base class and implements
    the conversion logic for reciprocal translocation.
    """

    def convert(self, pair):
        event1, event2 = pair
        # Check if this pair satisfies the criteria for reciprocal translocation
        if is_SpecialNoMateDiffBndPairReciprocalTranslocation(event1, event2):
            # Convert events
            self.convert_to_TRA(event1, event2)
            return [event1, event2]  # Return a list of transformed events
        else:
            return []  # If the pair doesn't satisfy the criteria, return an empty list

    def convert_to_TRA(self, event1, event2):
        # Convert a pair of events to reciprocal translocation
        # Modify event1
        event1.info["SVTYPE"] = "TRA"
        chrom_alt, pos_alt = get_alt_chrom_pos(event1.alt)
        event1.info["END"] = pos_alt
        event1.info["CHR2"] = chrom_alt
        event1.info["RTID"] = event2.id
        event1.info["SVLEN"] = 0
        # Modify event2
        event2.info["SVTYPE"] = "TRA"
        chrom_alt, pos_alt = get_alt_chrom_pos(event2.alt)
        event2.info["END"] = pos_alt
        event2.info["CHR2"] = chrom_alt
        event2.info["RTID"] = event1.id
        event2.info["SVLEN"] = 0


class SpecialNoMateDiffBNDPairIndependentToTRAConverter(Converter):
    """
    This class implements the conversion logic for independent translocation.
    """

    def convert(self, pair):
        event1, event2 = pair
        # Check if this pair satisfies the criteria for independent translocation
        if is_independent_special_bnd_event(event1, event2):
            # Convert events
            self.convert_to_TRA(event1, event2)
            return [event1, event2]  # Return a list of transformed events
        else:
            return []  # If the pair doesn't satisfy the criteria, return an empty list

    def convert_to_TRA(self, event1, event2):
        # Convert a pair of independent events to TRA
        # Modify event1
        event1.info["SVTYPE"] = "TRA"
        event1.info["SVLEN"] = 0
        chrom_alt, pos_alt = get_alt_chrom_pos(event1.alt)
        event1.info["END"] = pos_alt
        event1.info["CHR2"] = chrom_alt
        # Modify event2
        event2.info["SVTYPE"] = "TRA"
        event2.info["SVLEN"] = 0
        chrom_alt, pos_alt = get_alt_chrom_pos(event2.alt)
        event2.info["END"] = pos_alt
        event2.info["CHR2"] = chrom_alt


# ==============================
# The following strategies convert other_single_TRA to independent translocation.
# ==============================
class SingleTRAToTRAConverter(Converter):
    """
    This class inherits from the `Converter` base class and implements
    the conversion logic for single TRA events.
    """

    def convert(self, event):
        event.info["SVTYPE"] = "TRA"
        event.info["SVLEN"] = 0
        chrom_alt, pos_alt = get_alt_chrom_pos(event.alt)
        event.info["END"] = pos_alt
        event.info["CHR2"] = chrom_alt


# ==============================
# The following transformers receive different strategies and apply on the events and output the transformed events.
# ==============================


class EventTransformer:
    """Base class for transforming events."""

    def __init__(self, transform_strategies):
        self.transform_strategies = transform_strategies

    def apply_transforms(self, events):
        # The base implementation can be empty or provide a default behavior.
        pass

    def write_vcf(self, headers, transformed_events, output_file):
        # The base implementation can be empty or provide a default behavior.
        pass


class SameChrBNDTransformer(EventTransformer):  # The init input is list of strategies.
    """Class for transforming BND events on the same chromosome."""

    # It is initialized with a list of transform strategies,
    def apply_transforms(self, events):
        # Apply all transformation strategies to a list of events.
        for event in events:  # Try every converters for each event.
            for strategy in self.transform_strategies:
                strategy.convert(
                    event
                )  # Strategy is a converter instance, like BND_to_INV_Converter
        return events  # Returns: The transformed list of events.

    def write_vcf(self, headers, events, output_file):
        # Write the transformed events to a VCF file.
        with open(output_file, "w") as f:
            for header in headers:
                f.write(header + "\n")
            for event in events:
                f.write(str(event) + "\n")


class MatePairBNDTransformer(EventTransformer):
    """Class for transforming mate pair BND events."""

    def apply_transforms(self, mate_pairs):
        transformed_events = []
        for pair in mate_pairs:  # Pair is a tuple of mate events.
            for strategy in self.transform_strategies:
                transformed_events.extend(
                    strategy.convert(pair)
                )  # Strategy must be able to handle a tuple of events.
        return transformed_events


class SpecialNoMateDiffBNDPairTransformer(EventTransformer):
    """
    Class for transforming special_no_mate_diff_bnd_pair events.
    """

    def apply_transforms(self, event_pairs):
        transformed_events = []
        for pair in event_pairs:  # Pair is a tuple of two events.
            for strategy in self.transform_strategies:
                transformed_events.extend(
                    strategy.convert(pair)
                )  # Strategy must be able to handle a tuple of events.
        return transformed_events


class SingleTRATransformer(EventTransformer):
    """Class for transforming other_single_TRA."""

    def apply_transforms(self, events):
        for event in events:  # Apply the conversion for each event.
            for strategy in self.transform_strategies:
                strategy.convert(event)
        return events  # Return the original list of events, which have been modified in-place.


if __name__ == "__main__":
    main()

# Each software needs to calibrate the calculation method for SVLEN internally.
# After each software's conversion, redundancy must be removed internally.
# My standard format: Use <INV> for INV, use <DUP> for DUP, and for INS, DEL, it's best to keep the actual sequence.
# Now, only the non_bnd_events haven't been processed. They are actually more standard events like INS, INV, DEL, and DUP.
