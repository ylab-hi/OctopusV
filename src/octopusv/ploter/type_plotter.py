import collections
import logging

import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(level=logging.DEBUG)


class TypePlotter:
    def __init__(self, input_file):
        self.input_file = input_file
        # Define fixed SV type order and corresponding refined colors
        self.sv_order = ["TRA", "INV", "DUP", "INS", "DEL"]
        self.color_map = {
            "TRA": "#8e44ad",  # Muted purple
            "INV": "#c0392b",  # Deep crimson
            "DUP": "#d4ac0d",  # Soft gold
            "INS": "#27ae60",  # Forest green
            "DEL": "#2980b9",  # Ocean blue
        }
        self.data = self.parse_data()

    def parse_data(self):
        """Parse the statistics file and extract SV type information.

        Returns a dictionary with SV types as keys and (count, percentage) as values.
        """
        sv_types = {}
        parsing_types = False
        with open(self.input_file) as f:
            for line in f:
                if "SV Type Analysis" in line:
                    parsing_types = True
                    continue
                if parsing_types:
                    if line.strip() == "":
                        break
                    parts = line.strip().split("=")
                    if len(parts) == 2:
                        sv_type = parts[0].strip()
                        count_percentage = parts[1].strip().split()
                        if len(count_percentage) >= 2:
                            count = int(count_percentage[0])
                            percentage = float(count_percentage[1].strip("()%"))
                            sv_types[sv_type] = (count, percentage)

        logging.debug(f"Parsed SV type data: {sv_types}")
        return sv_types

    def _sort_data(self):
        """Sort data according to predefined order.
        Returns an OrderedDict with sorted SV types and their values.
        """
        ordered_data = collections.OrderedDict()
        # First add SV types in the predefined order
        for sv_type in self.sv_order:
            if sv_type in self.data:
                ordered_data[sv_type] = self.data[sv_type]

        # Add any additional SV types not in the predefined order at the end
        for sv_type, values in self.data.items():
            if sv_type not in self.sv_order:
                ordered_data[sv_type] = values

        return ordered_data

    def plot(self, output_prefix):
        """Create and save the SV type distribution plot.

        Args:
            output_prefix: Prefix for output files (.png and .svg will be appended)
        """
        if not self.data:
            logging.error("No data to plot")
            return

        # Sort data according to predefined order
        ordered_data = self._sort_data()

        # Set up the plot style
        plt.figure(figsize=(12, 8))
        sns.set_theme(style="white")

        types = list(ordered_data.keys())
        sizes = [count for count, _ in ordered_data.values()]

        # Get color list based on SV types
        colors = [self.color_map.get(sv_type, "#7f8c8d") for sv_type in types]  # Refined default gray

        # Create pie chart with refined aesthetics
        wedges, texts, autotexts = plt.pie(
            sizes,
            labels=types,
            autopct="%1.1f%%",
            startangle=90,
            colors=colors,
            wedgeprops={"width": 0.5, "edgecolor": "white", "linewidth": 2},
        )

        # Add center circle for donut chart with refined border
        centre_circle = plt.Circle((0, 0), 0.70, fc="white", ec="#f8f9fa")
        fig = plt.gcf()
        fig.gca().add_artist(centre_circle)

        # Ensure pie is drawn as a circle
        plt.axis("equal")

        # Refined title styling
        plt.title("SV Type Distribution", fontsize=18, pad=20, fontweight="bold", color="#2c3e50")

        # Create ordered legend labels with counts
        legend_labels = [f"{sv_type}: {ordered_data[sv_type][0]}" for sv_type in types]

        # Add legend with refined styling
        plt.legend(
            wedges,
            legend_labels,
            title="SV types",
            loc="center left",
            bbox_to_anchor=(1, 0, 0.5, 1),
            frameon=False,
            title_fontsize=12,
            fontsize=10,
        )

        # Adjust text properties for better readability
        for autotext in autotexts:
            autotext.set_fontsize(10)
            autotext.set_color("white")
            autotext.set_fontweight("bold")

        plt.tight_layout()
        # Save plots with high quality
        plt.savefig(f"{output_prefix}.png", dpi=300, bbox_inches="tight", facecolor="white")
        plt.savefig(f"{output_prefix}.svg", bbox_inches="tight", facecolor="white")
        plt.close()

        logging.info(f"Plot saved as {output_prefix}.png and {output_prefix}.svg")
