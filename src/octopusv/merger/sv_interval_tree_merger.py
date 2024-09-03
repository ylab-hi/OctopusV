from octopusv.merger.interval_tree import SVIntervalTree


class SVIntervalTreeMerger:
    """Base class for merging structural variant (SV) events from multiple sources."""

    def __init__(self):
        self.tree = SVIntervalTree()

    def load_events_into_tree(self, events, sv_type, chromosome):
        """Load SV events into the interval tree for a specific type and chromosome."""
        for event in events:
            start = event.pos
            end = event.end_position
            source_file = event.source_file
            self.tree.add_event(start, end, sv_type, chromosome, {source_file})

    def merge_by_type_and_chromosome(self):
        """Merge overlapping intervals within each type-specific and chromosome-specific tree."""
        self.tree.merge_overlaps()

    def write_result_by_type(self, output_file):
        """Write the results of merged SV events to a file."""
        with open(output_file, "w") as file:
            for sv_type, chrom_dict in self.tree.trees.items():
                for chrom, interval_tree in chrom_dict.items():
                    for interval in interval_tree:
                        file.write(f"{sv_type}, {chrom}: {interval.begin}-{interval.end}: {interval.data}\n")


# Example of structure before and after various operations

# After initialization:
"""
self.tree = {
    # No SV types or chromosomes added yet
}
"""

# After loading events:
"""
self.tree = {
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
            Interval(500, 600, {'fileD'}),
            Interval(450, 550, {'fileE'})
        ]),
        "chr8": IntervalTree([
            Interval(700, 800, {'fileF'})
        ])
    }
}
"""

# After merging overlaps:
"""
self.tree = {
    "INS": {
        "chr1": IntervalTree([
            Interval(100, 250, {'fileA', 'fileB'})  # Merged [100, 200) with [150, 250)
        ]),
        "chr5": IntervalTree([
            Interval(300, 400, {'fileC'})  # No overlaps
        ])
    },
    "DEL": {
        "chr5": IntervalTree([
            Interval(450, 600, {'fileD', 'fileE'})  # Merged [500, 600) with [450, 550)
        ]),
        "chr8": IntervalTree([
            Interval(700, 800, {'fileF'})  # No overlaps
        ])
    }
}
"""
