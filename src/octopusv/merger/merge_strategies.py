from octopusv.merger.interval_tree import SVIntervalTree


class MergeStrategy:
    """Base class for storing merge strategies."""

    def execute(self, trees):
        """Execute the merge strategy on a list of trees."""


class UnionStrategy(MergeStrategy):
    """Strategy to merge SV events by union."""

    def execute(self, trees):
        result_tree = SVIntervalTree()
        for tree in trees:
            for interval in tree.tree:
                result_tree.add_event(interval.begin, interval.end, *interval.data)
        result_tree.merge_overlaps()
        return result_tree


class IntersectionStrategy(MergeStrategy):
    """Strategy to merge SV events by intersection."""

    def execute(self, trees):
        if not trees:
            return None
        result_tree = trees[0]
        for tree in trees[1:]:
            result_tree = self._intersect_two(result_tree, tree)
        return result_tree

    def _intersect_two(self, tree1, tree2):
        """Helper function to perform intersection between two interval trees."""
        result_tree = SVIntervalTree()
        for interval1 in tree1.tree:
            overlaps = tree2.tree.overlap(interval1.begin, interval1.end)
            for interval2 in overlaps:
                start = max(interval1.begin, interval2.begin)
                end = min(interval1.end, interval2.end)
                if interval1.data[0] == interval2.data[0]:  # Ensure same type
                    result_tree.add_event(
                        start, end, interval1.data[0], interval1.data[1].intersection(interval2.data[1])
                    )
        return result_tree


class SpecificMergeStrategy(MergeStrategy):
    """Strategy to find common SV events across specified files."""

    def __init__(self, specific_files):
        self.specific_files = specific_files

    def execute(self, trees):
        """Execute the merge strategy to find common SV events in specified files."""
        common_tree = trees[self.specific_files[0]]
        for filename in self.specific_files[1:]:
            common_tree = self._intersect_two(common_tree, trees[filename])
        return common_tree

    def _intersect_two(self, tree1, tree2):
        """Intersect two interval trees and return the result."""
        result_tree = SVIntervalTree()
        for interval1 in tree1.tree:
            overlaps = tree2.tree.overlap(interval1.begin, interval1.end)
            for interval2 in overlaps:
                start = max(interval1.begin, interval2.begin)
                end = min(interval1.end, interval2.end)
                if interval1.data[0] == interval2.data[0]:  # Ensure same type
                    result_tree.add_event(start, end, interval1.data[0], interval1.data[1].union(interval2.data[1]))
        return result_tree
