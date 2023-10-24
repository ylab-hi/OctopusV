import re
import sys

from octopusv.sv import SVEvent


def is_same_chr_bnd(event):
    """Check if the POS and ALT of an event are on the same chromosome."""
    if event.is_BND():
        split_result = re.split(r"[\[\]:]", event.alt)
        if len(split_result) != 4:
            print(
                f"Unexpected ALT format, it should be something like N]chr10:69650962]: {split_result}",
            )
        else:
            chrom_alt, _ = split_result[1:3]
            return event.chrom == chrom_alt

    return False  # For non-BND, we won't categorize them as same_chr_bnd or diff_chr_bnd events


def check_vcf_format(vcf_file_path):
    """Check the format of a VCF file.

    Raise an error and exit if the format is incorrect.
    """
    with open(vcf_file_path) as f:
        lines = f.readlines()

    # Check if there is at least one header line
    if not any(line.startswith("#") for line in lines):
        print(
            "ERROR: Invalid VCF format. The file should contain at least one header line starting with '#'.",
        )
        sys.exit(1)

    for line in lines:
        if line.startswith("#"):
            continue  # Skip header lines

        # Check for space in lines
        if " " in line:
            print(
                "ERROR: Invalid VCF format. Non-header lines should not contain spaces.",
            )
            sys.exit(1)

        fields = line.split("\t")

        # Check the number of columns
        if len(fields) < 10:
            print(
                f"ERROR: Invalid VCF format. Expected at least 10 fields, but got {len(fields)}",
            )
            sys.exit(1)

        # Check that the position is a number
        try:
            int(fields[1])
        except ValueError:
            print(
                f"ERROR: Invalid VCF format. Position (field 2) should be a number, but got {fields[1]}",
            )
            sys.exit(1)

        # Check that the quality score is a number or '.'
        if fields[5] != ".":
            try:
                float(fields[5])
            except ValueError:
                print(
                    f"ERROR: Invalid VCF format. Quality score (field 6) should be a number or '.', but got {fields[5]}",
                )
                sys.exit(1)


def parse_vcf(vcf_file_path):
    """Parse VCF file into lists of SVEvent objects based on their type (same chromosome BND, different chromosome BND, non-BND) and return headers."""
    check_vcf_format(vcf_file_path)  # Check the format first
    same_chr_bnd_events = []
    diff_chr_bnd_events = []
    non_bnd_events = []
    headers = []

    with open(vcf_file_path) as f:
        for line in f:
            if line.startswith("#"):
                headers.append(line.strip())
                continue  # Skip header lines

            fields = line.strip().split("\t")
            event = SVEvent(*fields)  # Unpack fields and send to SVEvent class

            if event.is_BND():
                if is_same_chr_bnd(event):  # Check if the event is same chromosome
                    same_chr_bnd_events.append(event)
                else:
                    diff_chr_bnd_events.append(event)  # Different chromosome
            else:
                non_bnd_events.append(event)  # Non-BND events

    return headers, same_chr_bnd_events, diff_chr_bnd_events, non_bnd_events
