from octopusv.merger.interval_tree import SVIntervalTree

class SVMerger:
    """Base class for merging structural variant (SV) events from multiple sources."""
    def __init__(self):
        self.tree = SVIntervalTree()

    def load_events_into_tree(self, events):
        """Load SV events into the interval tree."""
        for event in events:
            start = event.pos
            end = event.end_position
            sv_type = event.sv_type
            source_file = event.source_file
            self.tree.add_event(start, end, sv_type, {source_file})

    def merge_by_type(self, sv_type, strategy):
        """Apply a merging strategy to the interval tree of a specific SV type."""
        specific_tree = self.tree.get_tree_by_type(sv_type)
        return strategy.execute([specific_tree])

    def write_result_by_type(self, output_file, sv_type):
        """Write the results of merged SV events of a specific type to a file."""
        with open(output_file, "w") as file:
            for interval in self.tree.get_tree_by_type(sv_type):
                file.write(f"{interval.begin}-{interval.end}: {interval.data}\n")
