# from elatool_arcpro import SingleDEM
from datetime import datetime
from enum import IntEnum

import arcpy


class Toolbox:
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "ELA Toolbox"
        self.alias = "ela_toolbox"
        self.description = "Contains tools for calculating ELAs on glacier surfaces and producing a variety of outputs."

        # List of tool classes associated with this toolbox
        self.tools = [SingleDEM]


class SingleDEM:

    class ParamIndex(IntEnum):
        DEM = 0
        METHOD = 1
        # OUT_CONTOURS = 2

    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = f"Single DEM TRIAL {datetime.now():%M:%S}"
        self.description = "Something descriptive here."

    def getParameterInfo(self) -> list[arcpy.Parameter]:
        """Define the tool parameters."""

        P = SingleDEM.ParamIndex
        params: list[arcpy.Parameter] = [arcpy.Parameter()]*len(P)

        params[P.DEM] = arcpy.Parameter(
            name="in_dem",
            displayName="Glacier surface DEM",
            datatype="DERasterDataset",
            multiValue=False,
            parameterType="Required",
            direction="Input")

        params[P.METHOD] = arcpy.Parameter(
            name="methods",
            displayName="ELA Calculation Methods",
            datatype="GPValueTable",
            # multiValue=True,
            parameterType="Required",
            direction="Input")
        params[P.METHOD].columns = [["GPString", "Method"],
                                    ["GPString", "Ratios"], ["GPString", "Intervals"]]
        params[P.METHOD].filters[0].type = "ValueList"
        params[P.METHOD].filters[0].list = ["AAR", "AABR"]

        return params

    def isLicensed(self):
        """Set whether the tool is licensed to execute."""
        return True

    def updateParameters(self, params: list[arcpy.Parameter]):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        P = SingleDEM.ParamIndex

        default_ratios = {
            "AAR": [0.50, 0.55, 0.58, 0.60, 0.65, 0.67, 0.70],
            "AABR": [1.0, 1.1, 1.56, 1.75, 2.0, 2.09, 2.47, 2.5, 3.0],
        }
        default_intervals = [100]

        if params[P.METHOD].value:
            value = params[P.METHOD].value
            for i, (method, ratios, intervals) in enumerate(value):
                if not ratios or len(str(ratios).strip()) == 0:
                    ratios = " ".join(str(m) for m in default_ratios[method])
                if not intervals or len(str(intervals).strip()) == 0:
                    intervals = " ".join(str(ival) for ival in default_intervals)
                value[i] = [method, ratios, intervals]
            params[P.METHOD].value = value

    def updateMessages(self, parameters: list[arcpy.Parameter]):
        """Modify the messages created by internal validation for each tool
        parameter. This method is called after internal validation."""
        pass

    def execute(self, params: list[arcpy.Parameter], messages):
        """The source code of the tool."""

        P = SingleDEM.ParamIndex

        arcpy.AddMessage("hello world")
        arcpy.AddMessage(params[P.DEM].value)
        arcpy.AddMessage(params[P.METHOD].value)

        return

    def postExecute(self, parameters: list[arcpy.Parameter]):
        """This method takes place after outputs are processed and
        added to the display."""
        pass
