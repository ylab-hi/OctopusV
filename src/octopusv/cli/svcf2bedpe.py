from pathlib import Path
import typer
from typing import Optional
from octopusv.utils.svcf_parser import SVCFFileEventCreator
from octopusv.formatter.svcf_to_bedpe_converter import SVCFtoBEDPEConverter

def svcf2bedpe(
    input_file: Path = typer.Option(..., "--input-file", "-i", exists=True, help="Input SVCF file to convert."),
    output_file: Path = typer.Option(..., "--output-file", "-o", help="Output BEDPE file."),
):
    """
    Convert SVCF file to BEDPE format.
    """
    # Parse SVCF file
    sv_event_creator = SVCFFileEventCreator(str(input_file))
    sv_event_creator.parse()

    # Convert to BEDPE
    converter = SVCFtoBEDPEConverter(sv_event_creator.events)
    bedpe_content = converter.convert()

    # Write to output file
    with open(output_file, 'w') as f:
        f.write(bedpe_content)

    typer.echo(f"Converted SVCF to BEDPE. Output written to {output_file}")

if __name__ == "__main__":
    typer.run(svcf2bedpe)