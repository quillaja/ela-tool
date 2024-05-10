from abc import ABC, abstractmethod
from typing import NamedTuple

import numpy as np


class HistData(NamedTuple):
    """Data returned by the Geoprocessor.histogram() function."""
    values: np.ndarray
    counts: np.ndarray
    bins: np.ndarray


class Geoprocessor(ABC):
    """
    This base class defines the interface for geoprocessing "back ends"
    required by the ELA tool.
    """

    # @abstractmethod
    # def surface_area(self, dem: str, side: Literal["ABOVE", "BELOW"], elevation: float) -> float:
    #     raise NotImplementedError

    @abstractmethod
    def message(self, msg: str) -> None:
        """
        Show a message via the Geoprocessor. For cli Geoprocessors, this
        may just wrap `print()`.

        :param msg: A message to be show.
        """
        raise NotImplementedError

    @abstractmethod
    def surface_area(self, dem: str, elevation: float) -> float:
        """
        Calculates the 3D surface area of the terrain in dem *below* the
        specified elevation level.

        :param dem: A path to a DEM of the surface.
        :param elevation: The elevation below which to find surface area.
        :return: The 3D surface area.
        """
        raise NotImplementedError

    @abstractmethod
    def min_max_elevations(self, dem: str) -> tuple[float, float]:
        """
        Get the minimum and maximum elevation values from dem.

        :param dem: A path to a DEM.
        :return: The min and max values as a tuple.
        """
        raise NotImplementedError

    @abstractmethod
    def array(self, dem: str) -> np.ndarray:
        """
        Get a numpy array of the values in dem. May be multidimensional.

        :param dem: A path to a DEM.
        :return: A numpy array of the elevation values.
        """
        raise NotImplementedError

    @abstractmethod
    def histogram(self, dem: str) -> HistData:
        """
        Get histogram data of dem. The returned values, counts, and bins
        are created by getting the DEM as an array using `array()`, flattening,
        removing NODATA values, creating 1 unit width bins from floor(min) to
        ceil(max). This is given to `np.histogram()` to generate counts and
        new bins.

        :param dem: A path to a DEM.
        :return: A NamedTuple containing the values, counts, and bins.
        """
        raise NotImplementedError
