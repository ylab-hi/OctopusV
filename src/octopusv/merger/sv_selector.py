def get_caller_score(caller_name, caller_priority_list):
    """
    Calculate a score for a structural variant caller based on its priority ranking.

    Args:
        caller_name (str): Name of the SV caller to evaluate
        caller_priority_list (list): Ordered list of caller names, where position indicates priority
                                (earlier position = higher priority)

    Returns:
        float: Score from 0-10, where:
            - Higher priority callers get higher scores
            - Score is normalized based on total number of known callers
            - Unknown callers are ranked last

    Example:
        If caller_priority_list = ['cuteSV', 'sniffles', 'SVIM'] and caller_name = 'cuteSV':
        - cuteSV is first (rank 1)
        - total_callers = 4 (3 known + 1 for unknown)
        - score = (4 - 1 + 1) / 4 * 10 = 10
    """
    try:
        rank = caller_priority_list.index(caller_name) + 1
    except ValueError:
        # If caller not in priority list, rank it last
        rank = len(caller_priority_list) + 1
    total_callers = len(caller_priority_list) + 1
    score = (total_callers - rank + 1) / total_callers * 10
    return score


def select_representative_sv(sv_group, weights=None, max_support=30, max_qual=100):
    """
    Select a representative SV record from a group of overlapping structural variants
    using a comprehensive scoring system.

    Args:
        sv_group (list): List of SVCFEvent objects representing overlapping structural variants
        weights (dict, optional): Dictionary specifying weights for different scoring components:
            - 'support': Weight for read support score (default: 0.5)
            - 'quality': Weight for quality score (default: 0.3)
            - 'caller': Weight for caller reliability score (default: 0.2)
        max_support (int, optional): Maximum read support value for normalization (default: 30)
        max_qual (float, optional): Maximum quality score for normalization (default: 100)

    Returns:
        SVCFEvent: The selected representative SV with merged source_file information

    Scoring System:
    1. Support Score (50% weight by default):
        - Based on number of supporting reads
        - Normalized to 0-10 range using max_support
        - Handles missing/invalid values gracefully

    2. Quality Score (30% weight by default):
        - Based on QUAL field in VCF
        - Normalized to 0-10 range using max_qual
        - Default value of 50 for missing/invalid scores

    3. Caller Score (20% weight by default):
        - Based on caller's reliability ranking
        - Predefined priority list: cuteSV > sniffles > SVIM > PBSV > ...
        - Unknown callers ranked last

    Tiebreaking Logic:
    If multiple SVs have the same total score, selection is based on:
    1. Higher read support
    2. If still tied, higher quality score

    Source File Handling:
    - Merges source_file information from all SVs in the group
    - Maintains sorted, unique list of sources
    """
    if weights is None:
        weights = {'support': 0.5, 'quality': 0.3, 'caller': 0.2}

    # Predefined priority list for SV callers based on empirical performance
    caller_priority_list = [
        'cuteSV',  # Best performance on long reads
        'sniffles',  # Well-established long read caller
        'SVIM',  # Good for PacBio data
        'PBSV',  # PacBio official caller
        'debreak',  # Newer caller with promising results
        'manta',  # Best short read caller
        'delly',  # Established short read caller
        'svaba',  # Specialized for somatic variants
        'svdss'  # Newer caller
    ]

    max_score = -1
    representative_sv = None
    all_source_files = set()  # Track all source files for merging

    for sv in sv_group:
        # Collect and update source file information
        source_files = sv.source_file.split(',')
        all_source_files.update(source_files)

        # Calculate support score with error handling
        support_value = sv.info.get('SUPPORT', '0')
        if support_value in [None, '.', '']:
            support = 0
        else:
            try:
                support = int(support_value)
            except ValueError:
                print(f"Warning: Invalid SUPPORT value '{support_value}' in SV {sv.sv_id}, setting support to 0.")
                support = 0
        support_score = min(support / max_support, 1.0) * 10

        # Calculate quality score with error handling
        qual_value = sv.quality
        if qual_value in [None, '.', '']:
            qual = 50  # Default quality score
        else:
            try:
                qual = float(qual_value)
            except ValueError:
                print(f"Warning: Invalid QUAL value '{qual_value}' in SV {sv.sv_id}, setting qual to 50.")
                qual = 50
        quality_score = min(qual / max_qual, 1.0) * 10

        # Calculate caller reliability score
        caller_name = sv.sample.get('SC', 'unknown')
        caller_score = get_caller_score(caller_name, caller_priority_list)

        # Calculate weighted total score
        total_score = (
                weights['support'] * support_score +
                weights['quality'] * quality_score +
                weights['caller'] * caller_score
        )

        # Update representative SV if current SV has higher score
        if total_score > max_score:
            max_score = total_score
            representative_sv = sv
        elif total_score == max_score:
            # Tiebreaker 1: Higher read support
            current_support = int(representative_sv.info.get('SUPPORT', '0')
                                if representative_sv.info.get('SUPPORT', '0') not in [None, '.', '']
                                else 0)
            if support > current_support:
                representative_sv = sv
            elif support == current_support:
                # Tiebreaker 2: Higher quality score
                current_qual = float(representative_sv.quality
                                    if representative_sv.quality not in [None, '.', '']
                                    else 50)
                if qual > current_qual:
                    representative_sv = sv

    # Merge source file information into representative SV
    representative_sv.source_file = ','.join(sorted(all_source_files))
    return representative_sv