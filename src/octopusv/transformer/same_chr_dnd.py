from .base import EventTransformer


class SameChrBNDTransformer(EventTransformer):  # The init input is list of strategies.
    """Class for transforming BND events on the same chromosome."""

    # It is initialized with a list of transform strategies,
    def apply_transforms(self, events):
        # Apply all transformation strategies to a list of events.
        for event in events:  # Try every converters for each event.
            for strategy in self.transform_strategies:
                strategy.convert(
                    event,
                )  # Strategy is a converter instance, like BND_to_INV_Converter
        return events  # Returns: The transformed list of events.

    def write_vcf(self, headers, events, output_file):
        # Write the transformed events to a VCF file.
        with open(output_file, "w") as f:
            for header in headers:
                f.write(header + "\n")
            for event in events:
                f.write(str(event) + "\n")
