from .base import SVMerger
from .merge_strategies import UnionStrategy, IntersectionStrategy

class SVCMerger(SVMerger):
    """
    Concrete implementation of the SVMerger to merge SV events from multiple files.
    """
    def merge(self, input_files, output_file, strategy):
        """
        Merge SV events using the specified strategy and write results to a file.
        """
        self.load_data(input_files)
        merged_tree = self.merge_by_type("any", strategy)  # Replace "any" with actual SV type
        self.write_result_by_type(output_file, "any")  # Replace "any" with actual SV type
