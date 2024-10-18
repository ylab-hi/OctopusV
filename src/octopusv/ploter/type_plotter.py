import matplotlib.pyplot as plt
import seaborn as sns
import logging

logging.basicConfig(level=logging.DEBUG)

class TypePlotter:
    def __init__(self, input_file):
        self.input_file = input_file
        self.data = self.parse_data()

    def parse_data(self):
        sv_types = {}
        parsing_types = False
        with open(self.input_file, 'r') as f:
            for line in f:
                if "SV Type Analysis" in line:
                    parsing_types = True
                    continue
                if parsing_types:
                    if line.strip() == "":
                        break
                    parts = line.strip().split('=')
                    if len(parts) == 2:
                        sv_type = parts[0].strip()
                        count_percentage = parts[1].strip().split()
                        if len(count_percentage) >= 2:
                            count = int(count_percentage[0])
                            percentage = float(count_percentage[1].strip('()%'))
                            sv_types[sv_type] = (count, percentage)
        logging.debug(f"Parsed SV type data: {sv_types}")
        return sv_types

    def plot(self, output_prefix):
        if not self.data:
            logging.error("No data to plot")
            return

        plt.figure(figsize=(12, 8))
        sns.set_style("whitegrid")
        colors = sns.color_palette("deep")

        types = list(self.data.keys())
        sizes = [count for count, _ in self.data.values()]
        percentages = [percentage for _, percentage in self.data.values()]

        # Create a pie chart
        wedges, texts, autotexts = plt.pie(sizes, labels=types, autopct='%1.1f%%', startangle=90, colors=colors,
            wedgeprops=dict(width=0.5, edgecolor='white'))

        # Add a circle at the center to create a donut chart
        centre_circle = plt.Circle((0,0), 0.70, fc='white')
        fig = plt.gcf()
        fig.gca().add_artist(centre_circle)

        # Equal aspect ratio ensures that pie is drawn as a circle
        plt.axis('equal')
        plt.title("SV Type Distribution", fontsize=18, pad=20)

        # Add legend with counts
        legend_labels = [f"{sv_type}: {count}" for sv_type, (count, _) in self.data.items()]
        plt.legend(wedges, legend_labels, title="SV types", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))

        # Adjust text size and color
        for autotext in autotexts:
            autotext.set_fontsize(10)
            autotext.set_color('white')

        plt.tight_layout()
        plt.savefig(f"{output_prefix}.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{output_prefix}.svg", bbox_inches='tight')
        plt.close()

        logging.info(f"Plot saved as {output_prefix}.png and {output_prefix}.svg")