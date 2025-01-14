import base64
import sys
import tempfile
from datetime import datetime
from pathlib import Path

from jinja2 import Template

from octopusv import __PACKAGE_NAME__
from octopusv.ploter.chromosome_plotter import ChromosomePlotter
from octopusv.ploter.size_plotter import SizePlotter
from octopusv.ploter.type_plotter import TypePlotter


def load_template() -> Path:
    """Load HTML template for report generation."""
    template_path = Path(sys.modules[__PACKAGE_NAME__].__file__)
    if template_path.parent is not None:
        template_path = template_path.parent / "report"
        return template_path / "template.html"
    raise FileNotFoundError("Could not locate template directory")


def image_to_base64(image_path):
    """Convert image to base64 string."""
    with Path.open(image_path, "rb") as image_file:
        try:
            return base64.b64encode(image_file.read()).decode('utf-8')
        except UnicodeDecodeError:
            # Try decoding with 'latin-1' if UTF-8 fails
            return base64.b64encode(image_file.read()).decode('latin-1')

class ReportGenerator:
    """Generate HTML report using Jinja2 template."""

    def __init__(self, template_path=None):
        """Initialize report generator with template path."""
        if template_path is None:
            # Use default template in same directory
            template_path = load_template()

        with Path.open(template_path) as f:
            self.template = Template(f.read())

    def plots(self, input_file, upset_plot=None):
        """Generate plots for the given input file."""
        chromosome_plotter = ChromosomePlotter(input_file)
        type_plotter = TypePlotter(input_file)
        size_plotter = SizePlotter(input_file)

        result = {}

        # generate three temperory files
        with tempfile.TemporaryDirectory() as tmpdirname:
            chromsome_plot = f"{tmpdirname}/chromosome_distribution.png"
            type_plot = f"{tmpdirname}/sv_types.png"
            size_plot = f"{tmpdirname}/sv_sizes.png"

            chromosome_plotter.plot(chromsome_plot)
            type_plotter.plot(type_plot)
            size_plotter.plot(size_plot)

            result['chromosome_coverage_plot'] = image_to_base64(chromsome_plot)
            result['size_distribution_plot'] = image_to_base64(size_plot)
            result['type_plot'] = image_to_base64(type_plot)

        if upset_plot:
            result['upset_plot'] = image_to_base64(upset_plot)

        return result


    def generate(self, input_file, output_path, sample_id):
        """Generate HTML report.

        Args:
            output_path: Path to save HTML report
            sample_id: Sample identifier
            summary_stats: Dict of summary statistics
            plots: Dict containing paths to generated plots:
                {
                    'sv_distribution_plot': path,
                    'chromosome_coverage_plot': path,
                    'size_distribution_plot': path,
                    'upset_plot': path,
                    'additional_plots': [
                        {'path': path, 'title': title, 'description': desc},
                        ...
                    ]
                }
        """
        plots = self.plots(input_file)

        # When preparing template data:
        template_data = {
            "generation_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "sample_id": sample_id,
            "chromosome_coverage_plot_base64": plots.get("chromosome_coverage_plot", ""),
            'sv_distribution_plot_base64': plots.get("type_plot", ""),
            "size_distribution_plot_base64": plots.get("size_distribution_plot", ""),
            'additional_plots': [{
                'title': plot.title,
                'description': plot.description,
                'base64_data': image_to_base64(plot.path)
            } for plot in plots]
        }

        # Render template
        html_content = self.template.render(**template_data)

        # Write output file
        with Path.open(output_path, "w") as f:
            f.write(html_content)
