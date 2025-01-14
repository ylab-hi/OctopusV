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

def load_logo() -> Path:
    """Load logo for report generation."""
    logo_path = Path(sys.modules[__PACKAGE_NAME__].__file__)
    if logo_path.parent is not None:
        logo_path = logo_path.parent / "report"
        return logo_path / "logo.png"
    raise FileNotFoundError("Could not locate logo directory")

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
            chromosome_plot_prefix = f"{tmpdirname}/chromosome_distribution"
            type_plot_prefix = f"{tmpdirname}/sv_types"
            size_plot_prefix = f"{tmpdirname}/sv_sizes"

            chromosome_plotter.plot(chromosome_plot_prefix, save_svg=False)
            type_plotter.plot(type_plot_prefix, save_svg=False)
            size_plotter.plot(size_plot_prefix, save_svg=False)

            chromosome_plot_path = f"{chromosome_plot_prefix}.png"
            type_plot_path = f"{type_plot_prefix}.png"
            size_plot_path = f"{size_plot_prefix}.png"

            if Path(chromosome_plot_path).exists():
                result['chromosome_coverage_plot'] = image_to_base64(chromosome_plot_path)

            if Path(size_plot_path).exists():
                result['size_distribution_plot'] = image_to_base64(size_plot_path)

            if Path(type_plot_path).exists():
                result['type_plot'] = image_to_base64(type_plot_path)

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
            "logo_path": load_logo(),
            "chromosome_coverage_plot_base64": plots.get("chromosome_coverage_plot", ""),
            'sv_distribution_plot_base64': plots.get("type_plot", ""),
            "size_distribution_plot_base64": plots.get("size_distribution_plot", ""),
        }

        # Render template
        html_content = self.template.render(**template_data)

        # Write output file
        with Path.open(output_path, "w") as f:
            f.write(html_content)
