from collections import Counter


class GenotypeAnalyzer:
    def __init__(self, input_file):
        self.input_file = input_file

    def analyze(self):
        genotypes = Counter()

        with open(self.input_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                format_fields = fields[8].split(":")
                sample_fields = fields[9].split(":")
                gt_index = format_fields.index("GT")
                gt = sample_fields[gt_index]
                genotypes[gt] += 1

        total = sum(genotypes.values())
        result = "Genotype Distribution:\n"
        for gt, count in genotypes.items():
            percentage = (count / total) * 100
            result += f"{gt}: {count} ({percentage:.2f}%)\n"

        return result
