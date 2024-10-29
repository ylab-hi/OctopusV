def should_merge_tra(event1, event2, delta=50, min_overlap_ratio=0.5, strand_consistency=True):
    # 检查染色体配对一致性
    if set([event1.start_chrom, event1.end_chrom]) != set([event2.start_chrom, event2.end_chrom]):
        return False

    # 检查方向性一致性
    if strand_consistency and not _is_compatible_bnd_pattern(event1.bnd_pattern, event2.bnd_pattern):
        return False

    # 计算断点位置差异
    pos1_diff = abs(event1.start_pos - event2.start_pos)
    pos2_diff = abs(event1.end_pos - event2.end_pos)
    if pos1_diff > 2 * delta or pos2_diff > 2 * delta:
        return False

    # 计算重叠区间的相似性
    overlap1 = max(
        0,
        min(event1.start_pos + delta, event2.start_pos + delta)
        - max(event1.start_pos - delta, event2.start_pos - delta),
    )
    overlap2 = max(
        0, min(event1.end_pos + delta, event2.end_pos + delta) - max(event1.end_pos - delta, event2.end_pos - delta)
    )

    ratio1 = overlap1 / (2 * delta)
    ratio2 = overlap2 / (2 * delta)

    return ratio1 >= min_overlap_ratio and ratio2 >= min_overlap_ratio


def _is_compatible_bnd_pattern(pattern1, pattern2):
    if pattern1 is None or pattern2 is None:
        return True  # 如果任一模式缺失，我们假设它们兼容
    if pattern1 == "<TRA>" and pattern2 == "<TRA>":
        return True

    def classify_pattern(pattern):
        if pattern == "<TRA>":
            return 5  # Special case for <TRA>
        if pattern.startswith("]") and pattern.endswith("N"):
            return 1  # ]chr:pos]N
        if pattern.startswith("N[") and pattern.endswith("["):
            return 2  # N[chr:pos[
        if pattern.startswith("N]") and pattern.endswith("]"):
            return 3  # N]chr:pos]
        if pattern.startswith("[") and pattern.endswith("N"):
            return 4  # [chr:pos[N
        return 0  # Unknown pattern

    class1 = classify_pattern(pattern1)
    class2 = classify_pattern(pattern2)

    return class1 == class2 and class1 != 0
