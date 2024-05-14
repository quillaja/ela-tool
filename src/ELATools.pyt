import time
from datetime import datetime
from enum import IntEnum
from itertools import chain
from pathlib import Path
from typing import Callable, Iterable, NamedTuple

import arcpy
import ela
from elatool_arcpro.esri import EsriGeoprocessor, elas_to_feature_class


class Toolbox:
    """The actual ArcPro toolbox. This cannot be renamed."""

    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "ELA Toolbox"
        self.alias = "ela_toolbox"
        self.description = "Contains tools for calculating ELAs on glacier surfaces and producing a variety of outputs."

        # List of tool classes associated with this toolbox
        self.tools = [BatchFindELATool]


def is_raster(p: str) -> bool:
    """Determine if `p` is a raster by trying to make an arcpy.Raster with it."""
    try:
        arcpy.Raster(p)
        return True
    except RuntimeError:
        pass
    return False


def is_workspace(p: str) -> bool:
    """A workspace is basically just a directory."""
    return Path(p).is_dir()


def workspace_rasters(ws: str) -> list[str]:
    """Get a list of paths to all rasters within workspace `ws`."""
    orig_ws = arcpy.env.workspace
    arcpy.env.workspace = ws
    rasters = [str(Path(ws, raster)) for raster in arcpy.ListRasters()]
    arcpy.env.workspace = orig_ws
    return rasters


def gather_rasters(semicolon_list: str) -> list[str]:
    """
    Split a list of paths to rasters or workspaces. If the path is to a workspace,
    get all rasters within that workspace. Combine all raster paths into
    one big list.
    """
    paths = semicolon_list.split(";")
    single_rasters = filter(is_raster, paths)
    ws_rasters = chain.from_iterable(workspace_rasters(ws)
                                     for ws in filter(is_workspace, paths))
    return list(chain(single_rasters, ws_rasters))


def split_strip(raw: str, sep: str = " ") -> Iterable[str]:
    """
    Split `raw` at `sep`, then strip whitespace from each part. Return
    only those parts that still have content.
    """
    for part in raw.strip().split(sep):
        part = part.strip()
        if part:
            yield part


class BatchFindELATool:
    """
    This tool will find ELAs for all combinations of glacier surface DEMs,
    ELA methods (AAR, AABR), intervals, and ratios. All the resultant ELAs are
    written to a new feature class along with contour lines at the ELA altitudes.
    Optionally, the data can be written to a CSV as well.
    """

    class MethodDefault(NamedTuple):
        ratios: list[float]
        intervals: list[float]
        factory: Callable[[float, str], ela.ELAMethod]

    class ParamIndex(IntEnum):
        DEM = 0
        METHOD = 1
        OUT_CONTOURS = 2
        OUT_CRS = 3
        OUT_CSV = 4

    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = f"Find ELAs DEM TRIAL {datetime.now():%H:%M:%S}"
        self.description = "Something descriptive here."

        gp = EsriGeoprocessor()

        self.defaults = {
            "AAR 3D": BatchFindELATool.MethodDefault(
                ratios=[0.50, 0.55, 0.58, 0.60, 0.65, 0.67, 0.70],
                intervals=[],
                factory=lambda ival, dem: ela.AAR(dem, gp)),

            "AAR 2D": BatchFindELATool.MethodDefault(
                ratios=[0.50, 0.55, 0.58, 0.60, 0.65, 0.67, 0.70],
                intervals=[],
                factory=lambda ival, dem: ela.AAR2D(dem, gp)),

            "AABR 3D": BatchFindELATool.MethodDefault(
                ratios=[1.0, 1.1, 1.56, 1.75, 2.0, 2.09, 2.47, 2.5, 3.0],
                intervals=[100],
                factory=lambda ival, dem: ela.AABR(ival, dem, gp)),

            "AABR 2D": BatchFindELATool.MethodDefault(
                ratios=[1.0, 1.1, 1.56, 1.75, 2.0, 2.09, 2.47, 2.5, 3.0],
                intervals=[],
                factory=lambda ival, dem: ela.AABR2D(dem, gp))
        }

    def getParameterInfo(self) -> list[arcpy.Parameter]:
        """Define the tool parameters."""

        P = BatchFindELATool.ParamIndex
        params: list[arcpy.Parameter] = [arcpy.Parameter()]*len(P)

        params[P.DEM] = arcpy.Parameter(
            name="in_dems",
            displayName="Glacier surface DEMs",
            datatype=["DERasterDataset", "Workspace"],
            multiValue=True,
            parameterType="Required",
            direction="Input")

        params[P.METHOD] = arcpy.Parameter(
            name="methods",
            displayName="ELA Calculation Methods",
            datatype="GPValueTable",
            multiValue=True,
            parameterType="Required",
            direction="Input")
        params[P.METHOD].columns = [("GPString", "Method"),
                                    ("GPString", "Ratios"),
                                    ("GPString", "Intervals")]
        params[P.METHOD].filters[0].type = "ValueList"
        params[P.METHOD].filters[0].list = list(self.defaults.keys())

        params[P.OUT_CONTOURS] = arcpy.Parameter(
            name="out_contours",
            displayName="Output ELA Contours Feature Class",
            datatype="DEFeatureClass",
            multiValue=False,
            parameterType="Required",
            direction="Output")
        params[P.OUT_CONTOURS].value = str(Path(arcpy.env.workspace, "ela_contours"))

        params[P.OUT_CRS] = arcpy.Parameter(
            name="out_crs",
            displayName="Output Spatial Reference",
            datatype="GPSpatialReference",
            multiValue=False,
            parameterType="Required",
            direction="Input")
        params[P.OUT_CRS].value = arcpy.mp.ArcGISProject("CURRENT").activeMap.spatialReference

        params[P.OUT_CSV] = arcpy.Parameter(
            name="out_csv",
            displayName="Output CSV",
            datatype="DETable",
            multiValue=False,
            parameterType="Optional",
            direction="Output",
            category="Additional Outputs")

        return params

    def isLicensed(self):
        """Set whether the tool is licensed to execute."""
        return True

    def updateParameters(self, params: list[arcpy.Parameter]):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        P = BatchFindELATool.ParamIndex

        # populate methods with default ratios and intervals
        if params[P.METHOD].value:
            value = params[P.METHOD].value
            for i, (method, ratios, intervals) in enumerate(value):
                if not ratios or len(str(ratios).strip()) == 0:
                    ratios = " ".join(str(m) for m in self.defaults[method].ratios)
                if not intervals or len(str(intervals).strip()) == 0:
                    intervals = " ".join(str(ival) for ival in self.defaults[method].intervals)
                value[i] = [method, ratios, intervals]
            params[P.METHOD].value = value

    def updateMessages(self, parameters: list[arcpy.Parameter]):
        """Modify the messages created by internal validation for each tool
        parameter. This method is called after internal validation."""
        pass

    def execute(self, params: list[arcpy.Parameter], messages):
        """The source code of the tool."""

        P = BatchFindELATool.ParamIndex

        # get inputs in sensible formats (mostly text)
        dem_raw: str = params[P.DEM].valueAsText
        methods: list[tuple[str, str, str]] = params[P.METHOD].value
        out_contours: str = params[P.OUT_CONTOURS].valueAsText
        out_crs: int = int(params[P.OUT_CRS].value.factoryCode)
        out_csv: str = params[P.OUT_CSV].valueAsText

        arcpy.AddMessage(f"{dem_raw=}\n{methods=}\n{out_contours=}\n{out_crs=}\n{out_csv=}")
        all_dems = gather_rasters(dem_raw)
        arcpy.AddMessage(f"{all_dems=}")

        elas: list[ela.ELA] = []
        for dem in all_dems:
            arcpy.AddMessage(f"processing {dem} ...")
            for method, ratios_raw, intervals_raw in methods:
                ratios: list[float] = [float(r) for r in split_strip(ratios_raw)]
                intervals: list[float] = [float(ival) for ival in split_strip(intervals_raw)]
                if len(intervals) == 0:
                    intervals.append(0)
                for ival in intervals:
                    start = time.perf_counter()
                    ela_method: ela.ELAMethod = self.defaults[method].factory(ival, dem)
                    elas += ela_method.estimate_ela(*ratios)
                    took = time.perf_counter()-start
                    arcpy.AddMessage(f"  {method} with interval {ival} took {took:.3f}s")

        # make feature class with ela results
        elas_to_feature_class(elas=elas,
                              out_fc=out_contours,
                              crs=out_crs)

        # write ela results to CSV if a path is specified
        if out_csv:
            ela.to_csv(out_csv, elas)

    def postExecute(self, parameters: list[arcpy.Parameter]):
        """This method takes place after outputs are processed and
        added to the display."""
        pass
