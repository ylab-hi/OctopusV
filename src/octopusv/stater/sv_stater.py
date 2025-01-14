from pathlib import Path

from .chromosome_analyzer import ChromosomeAnalyzer
from .genotype_analyzer import GenotypeAnalyzer
from .qc_analyzer import QCAnalyzer
from .size_analyzer import SizeAnalyzer
from .type_analyzer import TypeAnalyzer


class SVStater:
    def __init__(self, input_file, min_size=50, max_size=None):
        self.input_file = input_file
        self.min_size = min_size
        self.max_size = max_size
        self.results = {}

    def analyze(self):
        self.results["type"] = TypeAnalyzer(self.input_file).analyze()
        self.results["size"] = SizeAnalyzer(self.input_file, self.min_size, self.max_size).analyze()
        self.results["chromosome"] = ChromosomeAnalyzer(self.input_file).analyze()
        self.results["qc"] = QCAnalyzer(self.input_file).analyze()
        self.results["genotype"] = GenotypeAnalyzer(self.input_file).analyze()

    def write_results(self, output_file):
        with Path.open(output_file, "w") as f:
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
        for line in result.split("\n"):
            if ":" in line:
                key, value = line.split(":", 1)
                formatted += f"{key.strip():30} = {value.strip()}\n"
            else:
                formatted += line + "\n"
        return formatted

    def export_html(self, output_stat):
        """Export the analysis results to an HTML file."""
        result = {}
        result["input_file"] = self.input_file
        result["output_file"] = output_stat

        result['sv_types'] = {}
        for line in self.results["type"].strip().split("\n")[1:]:
            key, value = line.split(":")
            value_1, value_2 = value.strip().split(" ")
            value_2 = float(value_2[1:-2])
            result['sv_types'][key.strip()] = (int(value_1), value_2)

        lines = self.results["size"].strip().split("\n")
        result["total_svs"] = int(lines[0].split(":")[1].strip())
        result["min_size"] = int(lines[1].split(":")[1].strip().split()[0])
        result["max_size"] = int(lines[2].split(":")[1].strip().split()[0])
        result["mean_size"] = float(lines[3].split(":")[1].strip().split()[0])
        result["median_size"] = float(lines[4].split(":")[1].strip().split()[0])
        result['std_dev'] = float(lines[5].split(":")[1].strip().split()[0])

        result['size_distribution'] = {}
        for line in lines[8:]:
            key, value = line.split(":")
            result["size_distribution"][key.strip()] = int(value.strip())

        lines = self.results['qc'].strip().split("\n")
        result['avg_qual'] = float(lines[2].split(":")[1].strip())

        result["filter_status"] = {}
        pass_info = lines[10].split(":")[1].strip().split()
        result["filter_status"]["PASS"] = (pass_info[0], float(pass_info[1][1:-2]))
        hom_ref_info = lines[11].split(":")[1].strip().split()
        result["filter_status"]["hom_ref"] = (hom_ref_info[0], float(hom_ref_info[1][1:-2]))
        not_fully_covered_info = lines[12].split(":")[1].strip().split()
        result["filter_status"]["not_fully_covered"] = (not_fully_covered_info[0], float(not_fully_covered_info[1][1:-2]))

        result['avg_read_support'] = float(lines[16].split(":")[1].strip())


        lines = self.results['genotype'].strip().split("\n")
        result["genotype_dist"] = {}
        for line in lines[1:]:
            key, value = line.split(":")
            result["genotype_dist"][key.strip()] = (int(value.strip().split()[0]), float(value.strip().split()[1][1:-2]))

        return result

