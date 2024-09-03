from intervaltree import IntervalTree

"""
After add_event()

{
    "INS": {
        "chr1": IntervalTree([
            Interval(100, 200, {'fileA'}),
            Interval(150, 250, {'fileB'})
        ]),
        "chr5": IntervalTree([
            Interval(300, 400, {'fileC'})
        ])
    },
    "DEL": {
        "chr5": IntervalTree([
            Interval(500, 600, {'fileD'})
        ]),
        "chr8": IntervalTree([
            Interval(450, 550, {'fileE'})
        ])
    }
}
-----
After merge_overlaps()
{
    "INS": {
        "chr1": IntervalTree([
            Interval(100, 250, {'fileA', 'fileB'})  # Merge [100, 200)and[150, 250)
        ]),
        "chr5": IntervalTree([
            Interval(300, 400, {'fileC'})
        ])
    },
    "DEL": {
        "chr5": IntervalTree([
            Interval(500, 600, {'fileD'})
        ]),
        "chr8": IntervalTree([
            Interval(450, 550, {'fileE'})
        ])
    }
}


"""


class SVIntervalTree:
    """Class to manage SV events using interval trees, organized by SV type."""

    def __init__(self):
        self.trees = {}  # It's a dic, key is SVtype, value is another dic: key is chr, value are interval trees.

    def add_event(self, start, end, sv_type, chromosome, source):
        """Add an SV event to the appropriate type- and chromosome-specific interval tree."""
        if sv_type not in self.trees:
            self.trees[sv_type] = {}
        if chromosome not in self.trees[sv_type]:
            self.trees[sv_type][chromosome] = IntervalTree()
        self.trees[sv_type][chromosome].addi(start, end, source)

    def get_tree_by_type_and_chromosome(self, sv_type, chromosome):
        """Retrieve the interval tree for a specific SV type and chromosome."""
        return self.trees.get(sv_type, {}).get(chromosome, IntervalTree())

    def merge_overlaps(self):
        """Merge overlapping intervals within each type-specific and chromosome-specific tree."""
        for sv_type_trees in self.trees.values():
            for tree in sv_type_trees.values():
                tree.merge_overlaps(data_reducer=self._merge_source_files)

    def _merge_source_files(self, existing_sources, new_source):  # This function act on source file
        """Reducer function to combine sources of overlapping intervals."""
        return existing_sources.union(new_source)
