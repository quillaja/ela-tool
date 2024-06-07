from elatool_arcpro import (BatchFindELATool, CreateHistogramsTool,
                            MultiExtractRasterTool)


class Toolbox:
    """The actual ArcPro toolbox. This cannot be renamed."""

    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "ELA Toolbox"
        self.alias = "ela_toolbox"
        self.description = "Contains tools for calculating ELAs on glacier surfaces and producing a variety of outputs."

        # List of tool classes associated with this toolbox
        self.tools = [BatchFindELATool, MultiExtractRasterTool, CreateHistogramsTool]
