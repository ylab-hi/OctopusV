import sys

import typer

from src.octopusv import __version__

from .convert import convert

app = typer.Typer(
    epilog=f"{typer.style('Agent Octopus Code V helps you dive deep into the structural variations ocean!', fg=typer.colors.GREEN, bold=True)}",
    context_settings={"help_option_names": ["-h", "--help"]},
)


@app.command()
def convert_cmd():
    """Command to initiate convert functionality."""
    convert()

@app.callback()
def display_version_info():
    """Display the version information."""
    pass

display_version_info.__doc__ = f"""{typer.style("Version", fg=typer.colors.YELLOW, bold=True)}: {typer.style(f"{__version__}", fg=typer.colors.GREEN, bold=True)}"""


# For documentation purposes
if "sphinx" in sys.modules and __name__ != "__main__":
    # Create the typer click object to generate docs with sphinx-click
    typer_click_object = typer.main.get_command(app)


if __name__ == "__main__":
    app()
