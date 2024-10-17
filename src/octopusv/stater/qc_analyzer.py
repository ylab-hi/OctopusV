import statistics
from collections import defaultdict

class QCAnalyzer:
    def __init__(self, input_file):
        self.input_file = input_file

    def analyze(self):
        qual_scores = []
        filter_status = defaultdict(int)
        support_counts = []

        with open(self.input_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                qual = float(fields[5])
                filter_status[fields[6]] += 1

                info = dict(item.split('=') for item in fields[7].split(';'))
                support = int(info.get('SUPPORT', 0))

                qual_scores.append(qual)
                support_counts.append(support)

        result = "Quality Score (QUAL) Statistics:\n"
        result += f"Average QUAL: {statistics.mean(qual_scores):.2f}\n"
        result += f"Median QUAL: {statistics.median(qual_scores):.2f}\n"
        result += f"Min QUAL: {min(qual_scores):.2f}\n"
        result += f"Max QUAL: {max(qual_scores):.2f}\n"
        result += f"Q1 QUAL: {statistics.quantiles(qual_scores, n=4)[0]:.2f}\n"
        result += f"Q3 QUAL: {statistics.quantiles(qual_scores, n=4)[2]:.2f}\n\n"

        result += "Filter Status:\n"
        for status, count in filter_status.items():
            result += f"{status}: {count}\n"

        result += "\nRead Support Statistics:\n"
        result += f"Average Read Support: {statistics.mean(support_counts):.2f}\n"
        result += f"Median Read Support: {statistics.median(support_counts)}\n"
        result += f"Min Read Support: {min(support_counts)}\n"
        result += f"Max Read Support: {max(support_counts)}\n"
        result += f"Q1 Read Support: {statistics.quantiles(support_counts, n=4)[0]}\n"
        result += f"Q3 Read Support: {statistics.quantiles(support_counts, n=4)[2]}\n"

        return result