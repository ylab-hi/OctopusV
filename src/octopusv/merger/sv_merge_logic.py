def should_merge(event1, event2, max_distance=50, max_length_ratio=1.3, min_jaccard=0.7):
    """Determines whether two SV events should be merged.

    Enhanced version with multiple dimensions of analysis:
    1. Breakpoint consistency analysis
    2. Proportional analysis
    3. Position sensitivity analysis
    4. Local pattern analysis

    Args:
        event1: The first SV event
        event2: The second SV event
        max_distance: Maximum allowed distance between positions
        max_length_ratio: Maximum allowed ratio between lengths
        min_jaccard: Minimum required Jaccard index

    Returns:
        bool: True if events should be merged
    """
    # 1. Initial position checks
    start_diff = abs(event1.start_pos - event2.start_pos)
    end_diff = abs(event1.end_pos - event2.end_pos)

    if start_diff > max_distance or end_diff > max_distance:
        return False

    # 2. Length calculations
    length1 = event1.end_pos - event1.start_pos + 1
    length2 = event2.end_pos - event2.start_pos + 1
    min_length = min(length1, length2)
    max_length = max(length1, length2)

    # 3. Enhanced proportional analysis
    # Consider the proportion of difference relative to SV size
    start_diff_ratio = start_diff / min_length
    end_diff_ratio = end_diff / min_length

    # Reject if position difference is too large relative to SV size
    if start_diff_ratio > 0.2 or end_diff_ratio > 0.2:
        return False

    # 4. Position sensitivity analysis
    # Check if differences are proportional
    diff_ratio = abs(start_diff - end_diff) / max_distance
    if diff_ratio > 0.6:  # One end matches much better than the other
        return False

    # 5. Breakpoint consistency score
    bp_consistency = 1 - ((start_diff + end_diff) / (2 * max_distance))

    # 6. Enhanced length ratio analysis with size awareness
    length_ratio = max_length / min_length

    # Adjust ratio threshold based on SV size
    if min_length < 100:  # Small SVs
        adjusted_ratio = max_length_ratio * 0.9  # More strict
    elif min_length > 10000:  # Large SVs
        adjusted_ratio = max_length_ratio * 1.1  # More lenient
    else:  # Medium SVs
        # Linear scaling between 0.9 and 1.1
        scale_factor = (min_length - 100) / 9900  # 0 to 1
        adjusted_ratio = max_length_ratio * (0.9 + scale_factor * 0.2)

    if length_ratio > adjusted_ratio:
        return False

    # 7. Local pattern analysis
    # Check if the differences follow a consistent pattern
    position_pattern = abs(start_diff_ratio - end_diff_ratio)
    if position_pattern > 0.15:  # Differences should be somewhat consistent
        return False

    # 8. Enhanced overlap analysis
    overlap_start = max(event1.start_pos, event2.start_pos)
    overlap_end = min(event1.end_pos, event2.end_pos)

    if overlap_start <= overlap_end:
        overlap_length = overlap_end - overlap_start + 1
        union_start = min(event1.start_pos, event2.start_pos)
        union_end = max(event1.end_pos, event2.end_pos)
        union_length = union_end - union_start + 1

        # Basic Jaccard
        jaccard = overlap_length / union_length

        # Adjust required Jaccard based on consistency scores
        bp_factor = 0.9 + (bp_consistency * 0.2)  # 0.9 to 1.1
        size_factor = 0.9 + (min(length1, 10000) / 10000 * 0.2)  # 0.9 to 1.1

        required_jaccard = min_jaccard / (bp_factor * size_factor)

        # Additional check for small SVs
        if min_length < 100 and jaccard < required_jaccard * 1.1:
            return False

        return jaccard >= required_jaccard

    return False