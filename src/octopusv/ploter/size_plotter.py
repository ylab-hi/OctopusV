import logging

import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(level=logging.DEBUG)


class SizePlotter:
    def __init__(self, input_file):
        self.input_file = input_file
        self.data = self.parse_data()

    def parse_data(self):
        size_distribution = {}
        parsing_distribution = False
        with open(self.input_file) as f:
            for line in f:
                if "Size distribution" in line:
                    parsing_distribution = True
                    continue
                if parsing_distribution:
                    if line.strip() == "":
                        break
                    parts = line.strip().split("=")
                    if len(parts) == 2:
                        size_range = parts[0].strip()
                        count = int(parts[1].strip())
                        size_distribution[size_range] = count
        logging.debug(f"Parsed size distribution: {size_distribution}")
        return size_distribution

    def plot(self, output_prefix):
        if not self.data:
            logging.error("No data to plot")
            return

        plt.figure(figsize=(14, 8))
        sns.set_style("whitegrid")

        # Sort the size ranges
        size_order = ["0-50 bp", "51-100 bp", "101-500 bp", "501-1 kb", "1 kb-10 kb", ">10 kb"]
        sizes = [self.data.get(size, 0) for size in size_order]

        # Create bar plot
        bars = plt.bar(size_order, sizes, color=sns.color_palette("deep")[0], edgecolor="white", width=0.7)

        # Customize the plot
        plt.xlabel("SV Size Range", fontsize=14, labelpad=10)
        plt.ylabel("Count", fontsize=14, labelpad=10)
        plt.title("SV Size Distribution", fontsize=18, pad=20)
        plt.xticks(rotation=0, fontsize=12)
        plt.yticks(fontsize=12)

        # Add value labels on top of each bar
        for bar in bars:
            height = bar.get_height()
            plt.text(
                bar.get_x() + bar.get_width() / 2.0, height, f"{int(height)}", ha="center", va="bottom", fontsize=11
            )

        # Adjust layout and display
        plt.tight_layout()

        # Increase y-axis limit slightly to make room for value labels
        plt.ylim(0, max(sizes) * 1.1)

        # Add a light grid
        plt.grid(axis="y", linestyle="--", alpha=0.7)

        plt.savefig(f"{output_prefix}.png", dpi=300, bbox_inches="tight")
        plt.savefig(f"{output_prefix}.svg", bbox_inches="tight")
        plt.close()

        logging.info(f"Plot saved as {output_prefix}.png and {output_prefix}.svg")
