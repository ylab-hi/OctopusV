from .base import Converter


class NonBNDConverter(Converter):
    """A converter for non-BND events that doesn't change the event."""

    def convert(self, event):
        """No conversion needed for non-BND events for now will be revised soon."""
        pass
