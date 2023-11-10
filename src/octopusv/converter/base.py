import logging
import re

from natsort import natsorted

logging.basicConfig(level=logging.INFO)


class Converter:
    """This is an abstract base class for all converter classes.

    It provides a common interface for all converters.
    The `convert` method is a placeholder that needs to be overridden in each concrete converter class.
    """

    def convert(self, event):
        """This method should be overridden in a subclass."""
        raise NotImplementedError


def get_bnd_pattern(alt):
    """Get the pattern of BND: t[p[, t]p], ]p]t, [p[t."""
    if alt[0] in "ATCGN" and alt[1] == "[":
        return "t[p["
    if alt[0] in "ATCGN" and alt[1] == "]":
        return "t]p]"
    if alt[-1] in "ATCGN" and alt[-2] == "[":
        return "[p[t"
    if alt[-1] in "ATCGN" and alt[-2] == "]":
        return "]p]t"
    return None


def is_same_bnd_event(event1, event2) -> bool:
    """Define whether the BND represent the same TRA events, used by MatePairMergeToTRAConverter."""
    qualified_pairings = [{"]p]t", "t[p["}]
    event1_pattern = get_bnd_pattern(event1.alt)
    event2_pattern = get_bnd_pattern(event2.alt)

    return any({event1_pattern, event2_pattern} == pairing for pairing in qualified_pairings)


def is_independent_bnd_event(event1, event2):  # used by MatePairIndependentToTRAConverter
    """Determine if two mate_pair_bnd events are independent events."""
    qualified_pairings = [
        {"t]p]", "]p]t"},
        {"t]p]", "t]p]"},
        {"[p[t", "[p[t"},
        {"t[p[", "t]p]"},
        {"t[p[", "[p[t"},
        {"]p]t", "[p[t"},
    ]
    event1_pattern = get_bnd_pattern(event1.alt)
    event2_pattern = get_bnd_pattern(event2.alt)

    # I change back since tuple can not be compared with set.
    return any({event1_pattern, event2_pattern} == pairing for pairing in qualified_pairings)


def is_independent_special_bnd_event(
    event1,
    event2,
):  # Used by SpecialNoMateDiffBNDPairIndependentToTRAConverter
    """Determine if two special_no_mate_diff_bnd_pair events are independent events."""
    qualified_pairings = [
        {"t[p[", "t]p]"},  # Be sure to use set here.
        {"t[p[", "[p[t"},
        {"t]p]", "]p]t"},
        {"]p]t", "[p[t"},
    ]
    event1_pattern = get_bnd_pattern(event1.alt)
    event2_pattern = get_bnd_pattern(event2.alt)

    return any({event1_pattern, event2_pattern} == pairing for pairing in qualified_pairings)


def is_SpecialNoMateDiffBndPairReciprocalTranslocation(
    event1,
    event2,
):
    """Determine if SpecialNoMateDiffBndPairReciprocalTranslocation."""
    """used by SpecialNoMateDiffBNDPairReciprocalTranslocationToTRAConverter."""
    qualified_pairings = [{'t[p[', ']p]t'}, {'t]p]', '[p[t'}]
    # Extract the patterns from each event
    event1_pattern = get_bnd_pattern(event1.alt)
    event2_pattern = get_bnd_pattern(event2.alt)

    # Check if the patterns match one of the qualified pairings
    return any({event1_pattern, event2_pattern} == pairing for pairing in qualified_pairings)


def compare_chromosomes(event1, event2):  # Used by MatePairMergeToTRAConverter
    """Determine if the chromosome of event1 is smaller than the chromosome of event2 using natural sort order.

    Returns:
    A boolean value. True if the chromosome of event1 is smaller or equal to the chromosome of event2, False otherwise.
    """
    original_order = [event1.chrom, event2.chrom]

    # Get the order when the chromosomes are naturally sorted
    sorted_order = natsorted(original_order)

    # it means the chromosome of event1 is smaller or equal to the chromosome of event2
    return original_order == sorted_order


ALT_PARTS_EXPECTED_COUNT = 4  # Magic number for get_alt_chrom_pos()


def get_alt_chrom_pos(alt):
    """Extract chromosome and position from alt field in VCF file."""
    split_result = re.split(r"[\[\]:]", alt)
    if len(split_result) != ALT_PARTS_EXPECTED_COUNT:
        logging.info(
            f"Unexpected ALT format, it should be something like N]chr10:69650962]: {split_result}",
        )
        return None, None

    chrom_alt, pos_alt = split_result[1:3]
    return chrom_alt, int(pos_alt)  # Convert pos_alt to integer.
