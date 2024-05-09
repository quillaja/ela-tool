import math
import time
from abc import ABC, abstractmethod
from re import search
from typing import Iterable, NamedTuple

from networkx import maximal_matching

from .geoprocessing import Geoprocessor
from .ratio_calc import aar_direct, absolute_difference, br_from_slices
from .search import binary_search
from .slice import Slice, create_slices_threadpool


class ELA(NamedTuple):
    dem: str
    method: str
    interval: float
    ratio: float
    ela: float
    took: float


class ELAMethod(ABC):
    def __init__(self, dem: str, gp: Geoprocessor) -> None:
        super().__init__()
        self._gp = gp
        self.dem = dem

    @abstractmethod
    def estimate_ela(self, *ratios: float) -> list[ELA]:
        raise NotImplementedError()


class AAR(ELAMethod):
    def __init__(self, dem: str, gp: Geoprocessor) -> None:
        super().__init__(dem, gp)

    def estimate_ela(self, *ratios: float) -> list[ELA]:
        min_elev, max_elev = self._gp.min_max_elevations(self.dem)
        elas: list[ELA] = []
        for r in ratios:
            start = time.perf_counter()
            ela = binary_search(r, min_elev, max_elev,
                                aar_direct(self.dem),
                                absolute_difference(r))
            took = time.perf_counter() - start
            elas.append(ELA(
                dem=self.dem,
                method="AAR",
                interval=0,
                ratio=r,
                ela=ela,
                took=took))

        return elas


class AABR(ELAMethod):
    def __init__(self, interval: float, dem: str, gp: Geoprocessor) -> None:
        super().__init__(dem, gp)
        self.interval = interval
        self.slices: list[Slice] = []
        self.slices_took: float = 0

    def get_slices(self) -> list[Slice]:
        if not self.slices:
            start = time.perf_counter()
            self.slices = create_slices_threadpool(self.dem, self.interval)
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
                method="AABR",
                interval=self.interval,
                ratio=r,
                ela=ela,
                took=took+self.slices_took))

        return elas


class AABROriginal(ELAMethod):
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
