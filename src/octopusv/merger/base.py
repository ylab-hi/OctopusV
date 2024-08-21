class SVMerger:
    """
    Base class for merging structural variant (SV) events from multiple sources.
    """
    def __init__(self):
        self.tree = SVIntervalTree()

    def load_data(self, sv_files):
        """
        Load SV events from a list of files into the interval tree.
        Each file is expected to contain SV events in the format: type start end.
        """
        for filename in sv_files:
            with open(filename, 'r') as file:
                for line in file:
                    parts = line.strip().split()
                    start, end = map(int, parts[1:3])
                    sv_type = parts[0]
                    self.tree.add_event(start, end, sv_type, set([filename]))

    def merge_by_type(self, sv_type, strategy):
        """
        Apply a merging strategy to the interval tree of a specific SV type.
        """
        specific_tree = self.tree.get_tree_by_type(sv_type)
        return strategy.execute([specific_tree])

    def write_result_by_type(self, output_file, sv_type):
        """
        Write the results of merged SV events of a specific type to a file.
        """
        with open(output_file, 'w') as file:
            for interval in self.tree.get_tree_by_type(sv_type):
                file.write(f"{interval.begin}-{interval.end}: {interval.data}\n")
