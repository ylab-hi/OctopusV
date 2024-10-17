from .type_analyzer import TypeAnalyzer
from .size_analyzer import SizeAnalyzer
from .chromosome_analyzer import ChromosomeAnalyzer
from .qc_analyzer import QCAnalyzer
from .genotype_analyzer import GenotypeAnalyzer


class SVStater:
    def __init__(self, input_file, min_size=50, max_size=None):
        self.input_file = input_file
        self.min_size = min_size
        self.max_size = max_size
        self.results = {}

    def analyze(self):
        self.results['type'] = TypeAnalyzer(self.input_file).analyze()
        self.results['size'] = SizeAnalyzer(self.input_file, self.min_size, self.max_size).analyze()
        self.results['chromosome'] = ChromosomeAnalyzer(self.input_file).analyze()
        self.results['qc'] = QCAnalyzer(self.input_file).analyze()
        self.results['genotype'] = GenotypeAnalyzer(self.input_file).analyze()

    def write_results(self, output_file):
        with open(output_file, 'w') as f:
            f.write("OctopusV report\n")
            f.write("-" * 40 + "\n\n")

            f.write(">>>>>>> Input\n\n")
            f.write(f"input file = {self.input_file}\n")
            f.write(f"output file = {output_file}\n\n")

            for category, result in self.results.items():
                f.write(f">>>>>>> {category.capitalize()} Analysis\n\n")
                f.write(self.format_result(result))
                f.write("\n")

    def format_result(self, result):
        formatted = ""
        for line in result.split('\n'):
            if ':' in line:
                key, value = line.split(':', 1)
                formatted += f"{key.strip():30} = {value.strip()}\n"
            else:
                formatted += line + "\n"
        return formatted