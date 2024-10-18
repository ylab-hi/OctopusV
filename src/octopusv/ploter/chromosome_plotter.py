import matplotlib.pyplot as plt
import seaborn as sns
import logging

logging.basicConfig(level=logging.DEBUG)


class ChromosomePlotter:
    def __init__(self, input_file):
        self.input_file = input_file
        self.data = self.parse_data()

    def parse_data(self):
        chromosome_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
        chromosome_data = {chrom: 0 for chrom in chromosome_order}
        parsing_chromosomes = False
        with open(self.input_file, 'r') as f:
            for line in f:
                if "Chromosome Distribution" in line:
                    parsing_chromosomes = True
                    continue
                if parsing_chromosomes:
                    if line.strip() == "":
                        break
                    parts = line.strip().split('=')
                    if len(parts) == 2:
                        chrom = parts[0].strip()
                        if chrom in chromosome_data:
                            count = int(parts[1].split()[0])
                            chromosome_data[chrom] = count

        logging.debug(f"Parsed chromosome data: {chromosome_data}")
        return chromosome_data

    def plot(self, output_prefix):
        if not self.data:
            logging.error("No data to plot")
            return

        plt.figure(figsize=(20, 10))
        sns.set_style("whitegrid")
        colors = sns.color_palette("pastel")

        chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
        sv_counts = [self.data[chrom] for chrom in chromosomes]

        bars = plt.bar(chromosomes, sv_counts, color=colors[0])
        plt.xlabel("Chromosome", fontsize=14)
        plt.ylabel("Number of SVs", fontsize=14)
        plt.title("SV Distribution Across Chromosomes", fontsize=18)
        plt.xticks(rotation=0, ha='center', fontsize=12)
        plt.yticks(fontsize=12)

        # Add value labels on top of each bar
        for bar in bars:
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width() / 2., height, f'{int(height)}', ha='center', va='bottom', fontsize=10)

        # Adjust y-axis limit to add some space at the top
        plt.ylim(0, max(sv_counts) * 1.1)

        # Add a light grid
        plt.grid(axis='y', linestyle='--', alpha=0.7)

        plt.tight_layout()
        plt.savefig(f"{output_prefix}.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{output_prefix}.svg", bbox_inches='tight')
        plt.close()

        logging.info(f"Plot saved as {output_prefix}.png and {output_prefix}.svg")