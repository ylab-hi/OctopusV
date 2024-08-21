from intervaltree import IntervalTree, Interval

class SVIntervalTree:
    """
    Class to manage SV events using interval trees, organized by SV type.
    """
    def __init__(self):
        self.trees = {}

    def add_event(self, start, end, sv_type, source):
        """
        Add an SV event to the appropriate type-specific interval tree.
        """
        if sv_type not in self.trees:
            self.trees[sv_type] = IntervalTree()
        self.trees[sv_type].addi(start, end, source)

    def get_tree_by_type(self, sv_type):
        """
        Retrieve the interval tree for a specific SV type.
        """
        return self.trees.get(sv_type, IntervalTree())

    def merge_overlaps(self):
        """
        Merge overlapping intervals within each type-specific tree.
        """
        for tree in self.trees.values():
            tree.merge_overlaps(data_reducer=self._data_reducer)

    def _data_reducer(self, acc, item):
        """
        Reducer function to combine sources of overlapping intervals.
        """
        return acc.union(item)
