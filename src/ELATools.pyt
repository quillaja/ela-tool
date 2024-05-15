import time
from datetime import datetime
from enum import IntEnum
from itertools import chain
from pathlib import Path
from typing import Any, Callable, Iterable, NamedTuple, Optional

import arcpy
import arcpy.geoprocessing
import ela
from elatool_arcpro.esri import EsriGeoprocessor, elas_to_feature_class
from pytest import param


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
        self.label = f"Find ELAs {datetime.now():%H:%M:%S}"
        self.description = "The Find ELAs tool will find ELAs for all combinations of glacier surface DEMs,"\
            "ELA methods (AAR, AABR), intervals, and ratios. All the resultant ELAs are"\
            "written to a new feature class along with contour lines at the ELA altitudes."\
            "Optionally, the data can be written to a CSV as well."

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

    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = f"Extract Multiple Rasters {datetime.now():%H:%M:%S}"
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
            displayName="Python expression on 'name' to alter name field",
            datatype="GPString",
            multiValue=False,
            parameterType="Optional",
            direction="Input")
        params[P.EXPR].controlCLSID = "{E5456E51-0C41-4797-9EE4-5269820C6F0E}"  # multiline

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

        def is_layer(thing: Any) -> bool:
            try:
                thing.getSelectionSet()
                return True
            except AttributeError:
                pass
            return False

        def eval_name_expr(expr: str) -> Callable[[str], str]:
            return lambda name: str(eval(expr, {}, {"name": name}))

        def sql_select_where_from_selection(oids: set[str]) -> str:
            return f"OBJECTID IN ({','.join(str(id) for id in oids)})"

        arcpy.AddMessage(f"{dem=}\n{polygons=}\n{out_ws=}\n{name_field=}\n{name_expr=}")
        arcpy.AddMessage(polygons)

        clipping_fc: str = ""
        sql_clause: str = ""
        name_prep_func: Callable[[str], str] = None

        if is_layer(polygons):
            arcpy.AddMessage("i guess it's a Layer")
            clipping_fc = polygons.dataSource
            selection = polygons.getSelectionSet()
            sql_clause = sql_select_where_from_selection(selection) if selection is not None else ""
        else:
            arcpy.AddMessage("it's a string/'geoprocessing value object'")
            clipping_fc = str(polygons)

        if name_expr:
            name_prep_func = eval_name_expr(name_expr)

        multi_extract_raster(
            dem=dem,
            polygons=clipping_fc,
            output_location=out_ws,
            sql_select_where=sql_clause,
            name_field=name_field,
            name_prep=name_prep_func)

    def postExecute(self, parameters: list[arcpy.Parameter]):
        """This method takes place after outputs are processed and
        added to the display."""
        pass


def multi_extract_raster(dem: str, polygons: str, output_location: str,
                         sql_select_where: str = "",
                         name_field: Optional[str] = None,
                         name_prep: Optional[Callable[[str], str]] = None):

    fields = ["SHAPE@"]
    if name_field:
        fields.append(name_field)

    ext = ".tif" if Path(output_location).suffix.lower() != ".gdb" else ""

    with arcpy.da.SearchCursor(in_table=polygons,
                               field_names=fields,
                               where_clause=sql_select_where) as cursor:
        for i, row in enumerate(cursor):
            if name_field:
                if name_prep:
                    name = name_prep(str(row[1]))
                else:
                    name = str(row[1])
            else:
                name = f"{i:06d}"
            output = str(Path(output_location, name).with_suffix(ext))

            shape = row[0]
            clipped = arcpy.sa.ExtractByMask(in_raster=dem, in_mask_data=shape)
            clipped.save(output)
            arcpy.AddMessage(f"created {output}")
