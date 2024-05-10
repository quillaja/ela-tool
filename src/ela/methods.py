import math
import time
from abc import ABC, abstractmethod
from typing import NamedTuple

from .geoprocessing import Geoprocessor
from .ratio_calc import (aar_direct, aar_from_elev_hist, absolute_difference,
                         br_from_elev_hist, br_from_slices)
from .search import binary_search
from .slice import Slice, create_slices_threadpool


class ELA(NamedTuple):
    """
    ELA contains the result and parameters of an ELA calculation.
    """
    dem: str
    """Path to the DEM used for the ELA."""
    method: str
    """Method (eg AAR, AABR) used."""
    interval: float
    """The interval between vertical slices."""
    ratio: float
    """The ratio sought."""
    ela: float
    """The ELA found."""
    took: float
    """Processing time in seconds."""


class ELAMethod(ABC):
    """An ELAMethod defines the basic interface for by which to calculate ELAs."""

    def __init__(self, dem: str, gp: Geoprocessor) -> None:
        super().__init__()
        self._gp = gp
        self.dem = dem

    @abstractmethod
    def estimate_ela(self, *ratios: float) -> list[ELA]:
        """
        Perform the ELA calculation, producing one ELA result for each
        ratio in the `*ratios` argument.
        :param ratios: One or a list of ratios to find.
        :return: The ELAs found.
        """
        raise NotImplementedError()


class AAR(ELAMethod):
    """
    This ELAMethod uses binary search to quickly find the Accumulation Area
    Ratio method ELA.
    """

    def __init__(self, dem: str, gp: Geoprocessor) -> None:
        super().__init__(dem, gp)

    def estimate_ela(self, *ratios: float) -> list[ELA]:
        min_elev, max_elev = self._gp.min_max_elevations(self.dem)
        elas: list[ELA] = []
        for r in ratios:
            start = time.perf_counter()
            ela = binary_search(r, min_elev, max_elev,
                                aar_direct(self.dem, self._gp),
                                absolute_difference(r))
            took = time.perf_counter() - start
            elas.append(ELA(
                dem=self.dem,
                method="AAR_ben",
                interval=0,
                ratio=r,
                ela=ela,
                took=took))

        return elas


class AABR(ELAMethod):
    """
    This ELAMethod uses binary search with vertical slices of `interval` distance
    to find the Area Altitude Balance Ratio method ELA.
    """

    def __init__(self, interval: float, dem: str, gp: Geoprocessor) -> None:
        super().__init__(dem, gp)
        self.interval = interval
        self.slices: list[Slice] = []
        self.slices_took: float = 0

    def get_slices(self) -> list[Slice]:
        """Generate slices or retrieve cached slices."""
        if not self.slices:
            start = time.perf_counter()
            self.slices = create_slices_threadpool(self.dem, self.interval, self._gp)
            self.slices_took = time.perf_counter() - start
        return self.slices

    def estimate_ela(self, *ratios: float) -> list[ELA]:
        min_elev, max_elev = self._gp.min_max_elevations(self.dem)

        slices = self.get_slices()

        elas: list[ELA] = []
        for r in ratios:
            start = time.perf_counter()
            ela = binary_search(r, min_elev, max_elev,
                                br_from_slices(slices),
                                absolute_difference(r))
            took = time.perf_counter() - start
            elas.append(ELA(
                dem=self.dem,
                method="AABR_ben",
                interval=self.interval,
                ratio=r,
                ela=ela,
                took=took+self.slices_took))

        return elas


class AAR2D(ELAMethod):
    """
    This ELAMethod uses the DEM histogram to very quickly find ELAs using the
    Accumulation Area Ratio method based on 2D surface area.
    """

    def __init__(self, dem: str, gp: Geoprocessor) -> None:
        super().__init__(dem, gp)

    def estimate_ela(self, *ratios: float) -> list[ELA]:
        min_elev, max_elev = self._gp.min_max_elevations(self.dem)
        hist_data = self._gp.histogram(self.dem)

        elas: list[ELA] = []
        for r in ratios:
            start = time.perf_counter()
            ela = binary_search(r, min_elev, max_elev,
                                aar_from_elev_hist(hist_data.bins, hist_data.counts),
                                absolute_difference(r))
            took = time.perf_counter() - start
            elas.append(ELA(
                dem=self.dem,
                method="AAR_2D",
                interval=0,
                ratio=r,
                ela=ela,
                took=took))

        return elas


class AABR2D(ELAMethod):
    """
    This ELAMethod uses the DEM histogram to very quickly find ELAs using the
    Area Altitude Balance Ratio method based on 2D surface area.
    """

    def __init__(self, dem: str, gp: Geoprocessor) -> None:
        super().__init__(dem, gp)

    def estimate_ela(self, *ratios: float) -> list[ELA]:
        min_elev, max_elev = self._gp.min_max_elevations(self.dem)
        hist_data = self._gp.histogram(self.dem)

        elas: list[ELA] = []
        for r in ratios:
            start = time.perf_counter()
            ela = binary_search(r, min_elev, max_elev,
                                br_from_elev_hist(hist_data.bins, hist_data.counts),
                                absolute_difference(r))
            took = time.perf_counter() - start
            elas.append(ELA(
                dem=self.dem,
                method="AABR_2D",
                interval=0,
                ratio=r,
                ela=ela,
                took=took))

        return elas


class AABROriginal(ELAMethod):
    """
    This ELAMethod find the Area Altitude Balance Ratio using the 
    calculation presented by Pellitero, et al (2015) "A GIStool for automatic 
    calculation of glacier equilibrium-line altitudes". The code is mostly copied
    directly from Pellitero, but with a few modifications to make it function
    within this package's system.
    """

    def __init__(self, interval: float, dem: str, gp: Geoprocessor) -> None:
        super().__init__(dem, gp)
        self.interval = interval

    def estimate_ela(self, *ratios: float) -> list[ELA]:
        # 1. setup "slices"
        start = time.perf_counter()
        interval = int(self.interval)
        minalt, maxalt = self._gp.min_max_elevations(self.dem)
        maxalt = int(maxalt) + interval
        minalt = int(minalt) + interval

        # Create a list of altitudes
        list_altitudes: list[float] = []
        start_altitude = minalt+(interval/2)
        while start_altitude > minalt and start_altitude < maxalt:
            list_altitudes.append(start_altitude)
            start_altitude = start_altitude+interval

        # lists to be populated
        segundovalor: list[int] = []  # will hold surface area BELOW

        for plane in range(minalt, maxalt, interval):
            result = int(self._gp.surface_area(self.dem, plane))
            segundovalor.append(result)

        resta = [int(x)-int(y) for (x, y) in zip(segundovalor[1:], segundovalor[0:])]
        slices_took = time.perf_counter()-start

        # 2. find ELA from "slices" using ratio
        elas: list[ELA] = []
        for br in ratios:
            def aabr_calc(refinf: float) -> float:
                valores_multi = []
                valorAABR = [x*(y - refinf) for (x, y) in zip(resta, list_altitudes)]
                for valoracion in valorAABR:
                    if valoracion < 0:
                        valores_multi.append(valoracion*br)
                    else:
                        valores_multi.append(valoracion)
                valorAABRfinal = sum(valores_multi)
                return valorAABRfinal

            start = time.perf_counter()
            refinf = minalt-interval  # offset to make up for inital increment
            valorAABRfinal: float = math.inf
            while valorAABRfinal > 0:
                refinf += interval
                valorAABRfinal = aabr_calc(refinf)
            result = refinf - (interval/2)
            search_took = time.perf_counter() - start

            elas.append(ELA(
                dem=self.dem,
                method="AABR_original",
                interval=interval,
                ratio=br,
                ela=result,
                took=search_took+slices_took))

        return elas


class AAROriginal(ELAMethod):
    """
    This ELAMethod find the Accumulation Area Ratio using the 
    calculation presented by Pellitero, et al (2015) "A GIStool for automatic 
    calculation of glacier equilibrium-line altitudes". The code is mostly copied
    directly from Pellitero, but with a few modifications to make it function
    within this package's system.
    """

    def __init__(self, interval: float, dem: str, gp: Geoprocessor) -> None:
        super().__init__(dem, gp)
        self.interval = interval

    def estimate_ela(self, *ratios: float) -> list[ELA]:
        # 1. setup "slices"
        start = time.perf_counter()
        interval = int(self.interval)
        minalt, maxalt = self._gp.min_max_elevations(self.dem)
        maxalt = int(maxalt) + interval
        minalt = int(minalt) + interval

        # Create a list of altitudes
        primevalor: list[float] = []
        start_altitude = minalt  # NOTE: this didn't have + interval/2
        while start_altitude >= minalt and start_altitude < maxalt:
            primevalor.append(start_altitude)
            start_altitude = start_altitude+interval

        # lists to be populated
        segundovalor: list[int] = []  # will hold surface area ABOVE
        superf_total = self._gp.surface_area(self.dem, maxalt)  # total surface area
        for plane in range(minalt, maxalt, interval):
            result = int(superf_total - self._gp.surface_area(self.dem, plane))
            segundovalor.append(result)  # surface area ABOVE

        slices_took = time.perf_counter()-start

        elas: list[ELA] = []
        for ratio in ratios:
            start = time.perf_counter()
            # Create a list of the altitudes whose surface is less than ELA
            ELA_surface_area = superf_total*ratio
            superf_en_ELA: list[int] = []
            for values in segundovalor:
                if values <= ELA_surface_area:
                    superf_en_ELA.append(values)

            # Get the maximum surface value within the list
            ela = max(superf_en_ELA)

            # Get the altitude value that corresponds to the surface
            for altitudes, superficies in zip(primevalor, segundovalor):
                if superficies == ela:
                    elaok = altitudes + (interval/2)  # the answer

            search_took = time.perf_counter()-start

            elas.append(ELA(
                dem=self.dem,
                method="AAR_original",
                interval=interval,
                ratio=ratio,
                ela=elaok,
                took=slices_took+search_took))

        return elas
