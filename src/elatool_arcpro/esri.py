import math
import re
from functools import cache, reduce
from pathlib import Path
from typing import Callable, Iterable, Optional, Union

import arcpy
import ela
import numpy as np

# arcpy.CheckOutExtension("3D")  # do i need this?


class EsriGeoprocessor(ela.Geoprocessor):
    """A Geoprocessor() that interfaces with ArcGIS Pro via `arcpy`."""

    def message(self, msg: str) -> None:
        arcpy.AddMessage(msg)

    @cache
    def surface_area(self, dem: str, elevation: float) -> float:
        # https://pro.arcgis.com/en/pro-app/latest/tool-reference/3d-analyst/surface-volume.htm
        side = "BELOW"
        result = arcpy.SurfaceVolume_3d(dem, None, side, elevation).getMessage(1)
        # print(result)
        match = re.search(r"3D Area=\s*(?P<area>\d*[.,]?\d*)", result)
        area = match.group("area")
        return float(area)

    def min_max_elevations(self, dem: str) -> tuple[float, float]:
        # https://pro.arcgis.com/en/pro-app/latest/arcpy/classes/raster-object.htm
        r = arcpy.Raster(dem)
        return (r.minimum, r.maximum)

    # def cellsize(self, dem: str) -> tuple[float, float]:
    #     r = arcpy.Raster(dem)
        # #https://pro.arcgis.com/en/pro-app/latest/arcpy/classes/spatialreference.htm
        # crs: arcpy.SpatialReference = r.spatialReference
        # print(f"{crs.XYResolution = !s}")
        # print(f"{crs.abbreviation = !s}")
        # print(f"{crs.alias = !s}")
        # print(f"{crs.factoryCode = !s}")
        # print(f"{crs.name = !s}")
        # print(f"{crs.type = !s}")
        # print(f"{crs.linearUnitCode = !s}")
        # print(f"{crs.linearUnitName = !s}")
        # print(f"{crs.metersPerUnit = !s}")
        # print(f"{crs.projectionName = !s}")
        # print(f"{crs.datumName = !s}")
        # print(f"{crs.spheroidName = !s}")
        # if crs.VCS:
        #     print(f"{crs.VCS.datumName = !s}")
        #     print(f"{crs.VCS.factoryCode = !s}")
        #     print(f"{crs.VCS.linearUnitName = !s}")
        #     print(f"{crs.VCS.metersPerUnit = !s}")
        #     print(f"{crs.VCS.name = !s}")
        # else:
        #     print("no VCS")
        # return (r.meanCellWidth, r.meanCellHeight)

    # def avg_slope(self, dem: str) -> float:
    #     slope_raster = "memory\\"+Path(dem).stem+"_slope"
    #     arcpy.Slope_3d(dem, slope_raster, "DEGREE")
    #     slope = arcpy.Raster(slope_raster)
    #     arcpy.Delete_management(slope_raster)
    #     return slope.mean

    def array(self, dem: str) -> np.ndarray:
        # https://pro.arcgis.com/en/pro-app/latest/arcpy/functions/rastertonumpyarray-function.htm
        return arcpy.RasterToNumPyArray(dem, nodata_to_value=math.nan)

    def histogram(self, dem: str) -> ela.HistData:
        values: np.ndarray = self.array(dem)
        values = values.flatten()
        values = values[~np.isnan(values)]
        values.sort()
        low, high = math.floor(values.min()), math.ceil(values.max())
        bins = (x for x in range(low, high+1))
        counts, bins = np.histogram(values, bins=np.fromiter(bins, dtype=float))
        return ela.HistData(values, counts, bins)


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
                    arcpy.AddWarning(e)
