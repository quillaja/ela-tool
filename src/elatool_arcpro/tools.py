import time
from enum import IntEnum
from functools import reduce
from itertools import chain, groupby
from pathlib import Path
from typing import Any, Callable, Iterable, NamedTuple, Optional, Union

import arcpy
import ela
import numpy as np

from .esri import EsriGeoprocessor


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
                    took_nocache = sum(r.took for r in result)  # time without cache benefits
                    took_cache = time.perf_counter() - start  # real processing time
                    arcpy.AddMessage(
                        f"  {method} with interval {ival} and {len(ratios)} target ratios took {took_cache:.2f}s with cache, {took_nocache:.2f}s without cache.")

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


def elas_to_feature_class(elas: Iterable[ela.ELA], out_fc: str, crs: Union[int, arcpy.SpatialReference]) -> None:
    """
    Create or overwrite a new feature class at `out_fc` with the ELA
    information and contour lines contained within `elas`. The new feature class
    will be created with the spatial reference given in `crs`.

    :param elas: The ELA data to make into a feature class.
    :param out_fc: The feature class to create.
    :param crs: spatial reference for the new ELA feature class.
    """
    if arcpy.Exists(out_fc):
        arcpy.Delete_management(out_fc)

    ela_fields = [
        ["ela_id", "LONG"],
        ["ELA", "DOUBLE"],
        ["ratio", "DOUBLE"],
        ["interval", "DOUBLE"],
        ["method", "TEXT"],
        ["dem", "TEXT"]]
    ela_crs = crs if isinstance(
        crs, arcpy.SpatialReference) else arcpy.SpatialReference(crs)

    arcpy.CreateFeatureclass_management(
        out_path=str(Path(out_fc).parent),
        out_name=str(Path(out_fc).name),
        geometry_type="POLYLINE",
        spatial_reference=ela_crs)
    arcpy.AddFields_management(in_table=out_fc, field_description=ela_fields)

    # in order to add fields with ELA information, individual contours must be created
    # in memory, then their geometry can be written to a separate geodatabase along
    # with the field information. ContourList_3d() can write multiple contours at
    # once, but will not add a fields.
    temp_contour = r"memory\temp_contour"
    for id, r in enumerate(elas):
        arcpy.Delete_management(temp_contour)
        arcpy.ContourList_3d(r.dem, temp_contour, r.ela)
        with (arcpy.da.SearchCursor(temp_contour, ["SHAPE@"]) as source,
              arcpy.da.InsertCursor(out_fc, ["SHAPE@"]+[f[0] for f in ela_fields]) as dest):
            # make one multipart feature per ela
            geom: arcpy.Geometry = reduce(lambda a, b: a.union(b), [row[0] for row in source])
            dest.insertRow((geom, id, r.ela, r.ratio, r.interval, r.method, r.dem))
            # this part makes each contour part a new feature
            # for source_row in source:
            #     dest.insertRow((source_row[0], id, r.ela, r.ratio, r.interval, r.method, r.dem))


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


def multi_extract_raster(raster: str, polygons: str, output_location: str,
                         sql_select_where: Optional[str] = None,
                         name_field: Optional[str] = None,
                         name_prep: Optional[Callable[[str], str]] = None) -> None:
    """
    Uses the polygons to extract multiple rasters from a single raster. They
    will be saved to `output_location` as a arcpy `Raster` type if `output_location`
    is a geodatabase (.gdb) or as a TIFF (.tif) if `output_location` is a normal
    file folder.

    An sql "where" clause can be provided to select a subset of polygons. A field
    within polygons can be specified to use as the name of each new clipped raster,
    and the name can be transformed with the `name_prep` function. If no name_field
    is given, or if a row's value is NULL, the new rasters are named with
    sequential numbers from zero.

    :param raster: A path to a source raster that will be clipped.
    :param polygons: A path to a polygon feature class for clipping shapes.
    :param output_location: A path to a geodatabase or folder to put the new rasters.
    :param sql_select_where: A sql WHERE clause to get a subset of polygons.
    :param name_field: A field/attribute of polygons to use in naming the new rasters.
    :param name_prep: A function to transform name_field before saving the new rasters.
    """

    fields = ["SHAPE@"]
    if name_field:
        fields.append(name_field)

    ext = ".tif" if Path(output_location).suffix.lower() != ".gdb" else ""

    with arcpy.da.SearchCursor(in_table=polygons,
                               field_names=fields,
                               where_clause=sql_select_where) as cursor:
        for i, row in enumerate(cursor):
            if name_field and row[1]:
                if name_prep:
                    name = name_prep(str(row[1]))
                else:
                    name = str(row[1])
            else:
                raster_name = Path(raster).name
                name = f"{raster_name}_{i:06d}"
            output = str(Path(output_location, name).with_suffix(ext))

            shape = row[0]
            try:
                clipped = arcpy.sa.ExtractByMask(in_raster=raster, in_mask_data=shape)
                clipped.save(output)
                arcpy.AddMessage(f"created {output}")
            except arcpy.ExecuteError as e:
                if "010568" in str(e):  # 010568 is the error id as of 2024-05-15
                    arcpy.AddWarning(
                        "A polygon was skipped because it is outside the raster extents.")
                else:
                    raise e
            except RuntimeError as e:
                # RuntimeError:  ERROR 010240: on 2024-05-28 due to space in name of
                # raster saved to file geodatabase (.gdb)
                if "010240" in str(e):
                    arcpy.AddError(
                        f"Raster '{output}' contains a space or other character that is not allowed in the name of an item in a file geodatabase. Recommend adding 'name.replace(\" \", \"_\")' to the python expression field to customize output raster name.")
                else:
                    raise e


class CreateHistogramsTool:
    """
    The Create Histograms tool will produce surface area vs elevation
    histograms of glacier surfaces, and plot the ELAs on them.
    """

    class ParamIndex(IntEnum):
        ELAS = 0
        OUT_DIR = 1

    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = f"Create Histograms"
        self.description = "The Create Histograms tool will produce surface "\
            "area vs elevation histograpms of glacier surfaces."

    def getParameterInfo(self) -> list[arcpy.Parameter]:
        """Define the tool parameters."""

        P = CreateHistogramsTool.ParamIndex
        params: list[arcpy.Parameter] = [arcpy.Parameter()]*len(P)

        params[P.ELAS] = arcpy.Parameter(
            name="elas",
            displayName="ELAs to plot",
            datatype="GPFeatureLayer",
            multiValue=False,
            parameterType="Required",
            direction="Input")

        params[P.OUT_DIR] = arcpy.Parameter(
            name="out_dir",
            displayName="Output folder",
            datatype="DEFolder",
            multiValue=False,
            parameterType="Required",
            direction="Input")

        return params

    def isLicensed(self):
        """Set whether the tool is licensed to execute."""
        return True

    def updateParameters(self, params: list[arcpy.Parameter]):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        P = CreateHistogramsTool.ParamIndex

    def updateMessages(self, params: list[arcpy.Parameter]):
        """Modify the messages created by internal validation for each tool
        parameter. This method is called after internal validation."""

        P = CreateHistogramsTool.ParamIndex

    def execute(self, params: list[arcpy.Parameter], messages):
        """The source code of the tool."""

        P = CreateHistogramsTool.ParamIndex

        # get inputs in sensible formats (mostly text)
        elas: Any = params[P.ELAS].value  # arcpy.Value (?) or 'Layer'
        out_dir: str = params[P.OUT_DIR].valueAsText

        ela_fc: str = ""
        sql_clause: str = ""

        if is_layer(elas):
            # polygons is a "Layer" and may have an active selection
            ela_fc = elas.dataSource
            selection = elas.getSelectionSet()
            if selection:
                sql_clause = sql_select_where_from_selection(selection)
                arcpy.AddMessage("using only selected features")
        else:
            # polygons is a full path to a feature class
            ela_fc = str(elas)

        create_histograms(
            ela_fc=ela_fc,
            out_dir=out_dir,
            sql_clause=sql_clause)

    def postExecute(self, parameters: list[arcpy.Parameter]):
        """This method takes place after outputs are processed and
        added to the display."""
        pass


def create_histograms(ela_fc: str, out_dir: str, sql_clause: Optional[str] = None) -> None:
    gp = EsriGeoprocessor()

    ela_fields = ["ELA", "ratio", "interval", "method", "dem"]

    elas: list[ela.ELA] = []
    with arcpy.da.SearchCursor(ela_fc, ela_fields, where_clause=sql_clause) as cursor:
        for row in cursor:
            ela_elev, ratio, interval, method, dem = row
            elas.append(ela.ELA(dem, method, interval, ratio, ela_elev, 0))

    for (dem, method), grp in groupby(elas, lambda x: (x.dem, x.method)):
        grp = list(grp)

        out = Path(out_dir, f"{Path(dem).stem}_{method}.png")

        elevs = [e.ela for e in grp]
        ratios = [e.ratio for e in grp]
        ival = min(grp, key=lambda e: e.interval).interval
        ival = ival if ival > 0 else 50
        # arcpy.AddMessage(str(elevs))
        # arcpy.AddMessage(str(ratios))

        values: Union[np.ndarray, Iterable[ela.Slice]]
        if "2D" in method.upper():
            arcpy.AddMessage(f"getting 2D histogram of {dem}")
            ival = 1.0
            hist = gp.histogram(dem)
            values = hist.values

        elif "3D" in method.upper():
            arcpy.AddMessage(f"making slices for {dem} with interval {ival:.0f}")
            values = ela.make_slices(dem, ival, gp)

        else:
            arcpy.AddWarning(f"unknown ELA method for graphing: {method}")
            continue

        arcpy.AddMessage(f"saving {out}")
        ela.plot_histogram(
            data=values,
            ela=elevs,
            interval=ival,
            filename=str(out),
            title=f"Surface Area Dist. for {Path(dem).name} ({method}, interval={ival:.0f})")
