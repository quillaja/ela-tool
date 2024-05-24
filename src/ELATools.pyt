import time
from enum import IntEnum
from itertools import chain
from pathlib import Path
from typing import Any, Callable, Iterable, NamedTuple, Optional

import arcpy
import arcpy.geoprocessing
import ela
from elatool_arcpro.esri import (EsriGeoprocessor, elas_to_feature_class,
                                 multi_extract_raster)


class Toolbox:
    """The actual ArcPro toolbox. This cannot be renamed."""

    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "ELA Toolbox"
        self.alias = "ela_toolbox"
        self.description = "Contains tools for calculating ELAs on glacier surfaces and producing a variety of outputs."

        # List of tool classes associated with this toolbox
        self.tools = [BatchFindELATool, MultiExtractRasterTool]


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
    The Find ELAs tool will find ELAs for all combinations of glacier surface DEMs,
    ELA methods (AAR, AABR), intervals, and ratios. All the resultant ELAs are
    written to a new feature class along with contour lines at the ELA altitudes.
    Optionally, the data can be written to a CSV as well.
    """

    # see defaults dict in __init__
    class MethodDefault(NamedTuple):
        ratios: list[float]
        intervals: list[float]
        factory: Callable[[float, str], ela.ELAMethod]

    # indices for tool input parameters
    class ParamIndex(IntEnum):
        DEM = 0
        METHOD = 1
        OUT_CONTOURS = 2
        OUT_CRS = 3
        OUT_CSV = 4

    # this class variable is meant as a method "memo" so that in updateParameters()
    # we can check if the user changed the selected method or not. Apparently ArcGIS
    # creates new instances of this class all the damn time, so making this an instance
    # variable doesn't work (see link). The dict is reset in the execute() method.
    # https://gis.stackexchange.com/questions/156939/python-toolbox-tool-objects-do-not-remember-instance-variables-between-function
    _method_selection: dict[int, str] = {}

    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = f"Find ELAs"
        self.description = "The Find ELAs tool will find ELAs for all combinations of glacier surface DEMs,"\
            "ELA methods (AAR, AABR), intervals, and ratios. All the resultant ELAs are"\
            "written to a new feature class along with contour lines at the ELA altitudes."\
            "Optionally, the data can be written to a CSV as well."

        gp = EsriGeoprocessor()  # access to arcpy from within the ela package

        # this dictionary configures the geoprocessing tool's defaults for
        # ELA methods. It also provides a "factory" function to more easily
        # create new ELAMethod instances from parameters
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
                factory=lambda ival, dem: ela.AABR2D(dem, gp)),

            "AAR 3D Pellitero": BatchFindELATool.MethodDefault(
                ratios=[0.50, 0.55, 0.58, 0.60, 0.65, 0.67, 0.70],
                intervals=[100],
                factory=lambda ival, dem: ela.AAROriginal(ival, dem, gp)),

            "AABR 3D Pellitero": BatchFindELATool.MethodDefault(
                ratios=[1.0, 1.1, 1.56, 1.75, 2.0, 2.09, 2.47, 2.5, 3.0],
                intervals=[100],
                factory=lambda ival, dem: ela.AABROriginal(ival, dem, gp))
        }

    def getParameterInfo(self) -> list[arcpy.Parameter]:
        """Define the tool parameters."""

        P = BatchFindELATool.ParamIndex
        params: list[arcpy.Parameter] = [arcpy.Parameter()]*len(P)

        params[P.DEM] = arcpy.Parameter(
            name="in_dems",
            displayName="Glacier surface DEMs",
            datatype=["GPRasterLayer", "Workspace"],
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
        params[P.OUT_CRS].value = arcpy.SpatialReference(
            arcpy.mp.ArcGISProject("CURRENT").activeMap.spatialReference.factoryCode)

        params[P.OUT_CSV] = arcpy.Parameter(
            name="out_csv",
            displayName="Optional CSV output location",
            datatype="DETable",
            multiValue=False,
            parameterType="Optional",
            direction="Output")

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
                # check if current method is the method previously selected at this index,
                # and if not, reset the ratios and intervals and save the new method selection
                if i not in self._method_selection or method != self._method_selection[i]:
                    ratios = ""
                    intervals = ""
                    self._method_selection[i] = method
                # set default ratios and intervals for the selected method if there
                # is nothing already entered by the user
                if not ratios or len(str(ratios).strip()) == 0:
                    ratios = " ".join(str(m) for m in self.defaults[method].ratios)
                if not intervals or len(str(intervals).strip()) == 0:
                    intervals = " ".join(str(ival) for ival in self.defaults[method].intervals)
                # write the updated method parameters
                value[i] = [method, ratios, intervals]
            params[P.METHOD].value = value

    def updateMessages(self, parameters: list[arcpy.Parameter]):
        """Modify the messages created by internal validation for each tool
        parameter. This method is called after internal validation."""
        pass

    def execute(self, params: list[arcpy.Parameter], messages):
        """The source code of the tool."""

        self._method_selection.clear()  # reset class method memo when finally executing
        P = BatchFindELATool.ParamIndex

        # get inputs in sensible formats (mostly text)
        dem_raw: str = params[P.DEM].valueAsText
        methods: list[tuple[str, str, str]] = params[P.METHOD].value
        out_contours: str = params[P.OUT_CONTOURS].valueAsText
        out_crs: int = int(params[P.OUT_CRS].value.factoryCode)
        out_csv: str = params[P.OUT_CSV].valueAsText

        all_dems = gather_rasters(dem_raw)  # all possible input rasters to 1 master list

        # perform ELA calculation on each DEM, method, interval, and ratio.
        # ratios and intervals input will be cleaned/preprepared from strings
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
                    result = ela_method.estimate_ela(*ratios)
                    elas += result
                    took = (time.perf_counter()-start)/len(ratios)
                    arcpy.AddMessage(
                        f"  {method} with interval {ival} took {took:.3f}s per ratio")

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


def is_layer(thing: Any) -> bool:
    """Determine if the 'feature layer' is a Layer object."""
    try:
        thing.getSelectionSet()
        return True
    except AttributeError:
        pass
    return False


def eval_name_expr(expr: str) -> Callable[[str], str]:
    """Make a function that uses the python expression to transform a string."""
    return lambda name: str(eval(expr, {}, {"name": name}))


def sql_select_where_from_selection(oids: set[str]) -> str:
    """Prepare a sql statment with the set of object ids."""
    return f"OBJECTID IN ({','.join(str(id) for id in oids)})"


class MultiExtractRasterTool:
    """
    The Extract Multiple Rasters tool will extract multiple rasters from a single
    large raster using the clipping/masking polygons in an input feature class. The
    resultant rasters can be written to a geodatabase or folder.
    """

    class ParamIndex(IntEnum):
        DEM = 0
        POLYGONS = 1
        OUT_WS = 2
        NAME = 3
        EXPR = 4
        SQL = 5

    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = f"Extract Multiple Rasters"
        self.description = "The Extract Multiple Rasters tool will extract multiple rasters "\
            "from a single large raster using the clipping/masking polygons in an input "\
            "feature class. The resultant rasters can be written to a geodatabase or folder."

    def getParameterInfo(self) -> list[arcpy.Parameter]:
        """Define the tool parameters."""

        P = MultiExtractRasterTool.ParamIndex
        params: list[arcpy.Parameter] = [arcpy.Parameter()]*len(P)

        params[P.DEM] = arcpy.Parameter(
            name="in_dem",
            displayName="Source raster",
            datatype="GPRasterLayer",
            multiValue=False,
            parameterType="Required",
            direction="Input")

        params[P.POLYGONS] = arcpy.Parameter(
            name="polygons",
            displayName="Clipping polygons",
            datatype="GPFeatureLayer",
            multiValue=False,
            parameterType="Required",
            direction="Input")
        params[P.POLYGONS].filter.list = ["Polygon"]

        params[P.OUT_WS] = arcpy.Parameter(
            name="out_ws",
            displayName="Output workspace",
            datatype="DEWorkspace",
            multiValue=False,
            parameterType="Required",  # watch out setting this to derived
            direction="Input")
        params[P.OUT_WS].value = str(arcpy.env.workspace)

        params[P.NAME] = arcpy.Parameter(
            name="name_field",
            displayName="Name field in clipping polygons",
            datatype="Field",
            multiValue=False,
            parameterType="Optional",
            direction="Input")
        params[P.NAME].filter.list = ["Text"]
        params[P.NAME].parameterDependencies = [params[P.POLYGONS].name]

        params[P.EXPR] = arcpy.Parameter(
            name="name_expr",
            displayName="Python expression with 'name' to customize output raster name",
            datatype="GPString",
            multiValue=False,
            parameterType="Optional",
            direction="Input")
        params[P.EXPR].controlCLSID = "{E5456E51-0C41-4797-9EE4-5269820C6F0E}"  # multiline

        params[P.SQL] = arcpy.Parameter(
            name="sql_filter",
            displayName="Polygons filter SQL expression",
            datatype="GPSQLExpression",
            multiValue=False,
            parameterType="Optional",
            direction="Input")
        params[P.SQL].parameterDependencies = [params[P.POLYGONS].name]

        return params

    def isLicensed(self):
        """Set whether the tool is licensed to execute."""
        return True

    def updateParameters(self, params: list[arcpy.Parameter]):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        P = MultiExtractRasterTool.ParamIndex

    def updateMessages(self, params: list[arcpy.Parameter]):
        """Modify the messages created by internal validation for each tool
        parameter. This method is called after internal validation."""

        P = MultiExtractRasterTool.ParamIndex

        # check and set error if the python expression has content but
        # no name field was chosen
        if params[P.EXPR].valueAsText and (not params[P.NAME].valueAsText):
            params[P.NAME].setErrorMessage(
                "A name field must be chosen if a python expression is given.")

    def execute(self, params: list[arcpy.Parameter], messages):
        """The source code of the tool."""

        P = MultiExtractRasterTool.ParamIndex

        # get inputs in sensible formats (mostly text)
        dem: str = params[P.DEM].valueAsText
        polygons: Any = params[P.POLYGONS].value  # arcpy.Value (?) or 'Layer'
        out_ws: str = params[P.OUT_WS].valueAsText
        name_field: Optional[str] = params[P.NAME].valueAsText  # None if blank
        name_expr: Optional[str] = params[P.EXPR].valueAsText  # None if blank
        sql_clause: Optional[str] = params[P.SQL].valueAsText  # None if blank

        clipping_fc: str = ""
        name_prep_func: Optional[Callable[[str], str]] = None

        if is_layer(polygons):
            # polygons is a "Layer" and may have an active selection
            clipping_fc = polygons.dataSource
            selection = polygons.getSelectionSet()
            if selection:
                sql_selection = sql_select_where_from_selection(selection)
                sql_clause = f"{sql_clause} AND {sql_selection}" if sql_clause else sql_selection
        else:
            # polygons is a full path to a feature class
            clipping_fc = str(polygons)

        # if a python expression is provided, wrap the call to 'eval'
        # in a function that can be passed to multi_extract_raster()
        if name_expr:
            name_prep_func = eval_name_expr(name_expr)

        # do it
        multi_extract_raster(
            raster=dem,
            polygons=clipping_fc,
            output_location=out_ws,
            sql_select_where=sql_clause,
            name_field=name_field,
            name_prep=name_prep_func)

    def postExecute(self, parameters: list[arcpy.Parameter]):
        """This method takes place after outputs are processed and
        added to the display."""
        pass
