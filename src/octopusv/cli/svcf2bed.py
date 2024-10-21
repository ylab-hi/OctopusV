from pathlib import Path
import typer
from typing import Optional
from octopusv.utils.svcf_parser import SVCFFileEventCreator
from octopusv.formatter.svcf_to_bed_converter import SVCFtoBEDConverter

def svcf2bed(
    input_file: Path = typer.Option(..., "--input-file", "-i", exists=True, help="Input SVCF file to convert."),
    output_file: Path = typer.Option(..., "--output-file", "-o", help="Output BED file."),
):
    """
    Convert SVCF file to BED format.
    """
    # Parse SVCF file
    sv_event_creator = SVCFFileEventCreator(str(input_file))
    sv_event_creator.parse()

    # Convert to BED
    converter = SVCFtoBEDConverter(sv_event_creator.events)
    bed_content = converter.convert()

    # Write to output file
    with open(output_file, 'w') as f:
        f.write(bed_content)

    typer.echo(f"Converted SVCF to BED. Output written to {output_file}")

if __name__ == "__main__":
    typer.run(svcf2bed)