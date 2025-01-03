from datetime import datetime
from pathlib import Path

from jinja2 import Template


import base64

def image_to_base64(image_path):
    with open(image_path, "rb") as image_file:
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
            template_path = Path(__file__).parent / "template.html"

        with open(template_path) as f:
            self.template = Template(f.read())

    def generate(self, output_path, sample_id, summary_stats, plots):
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
                    'quality_metrics_plot': path,
                    'additional_plots': [
                        {'path': path, 'title': title, 'description': desc},
                        ...
                    ]
                }
        """
        # Prepare template variables

        # When preparing template data:
        template_data = {
            "generation_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "sample_id": sample_id,
            "summary_stats": summary_stats,
            'logo_base64': image_to_base64('path/to/logo.png'),
            'sv_distribution_plot_base64': image_to_base64('path/to/sv_plot.png'),
            "chromosome_coverage_plot_base64": plots.get("chromosome_coverage_plot", ""),
            "size_distribution_plot_base64": plots.get("size_distribution_plot", ""),
            "quality_metrics_plot_base64": plots.get("quality_metrics_plot", ""),
            'additional_plots': [{
                'title': plot.title,
                'description': plot.description,
                'base64_data': image_to_base64(plot.path)
            } for plot in plots]
        }

        # Render template
        html_content = self.template.render(**template_data)

        # Write output file
        with open(output_path, "w") as f:
            f.write(html_content)
