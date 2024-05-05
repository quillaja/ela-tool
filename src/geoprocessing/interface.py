from abc import ABC, abstractmethod
from typing import Literal, NamedTuple

import numpy as np


class HistData(NamedTuple):
    values: np.ndarray
    counts: np.ndarray
    bins: np.ndarray


class Geoprocessor(ABC):

    # @abstractmethod
    # def surface_area(self, dem: str, side: Literal["ABOVE", "BELOW"], elevation: float) -> float:
    #     raise NotImplementedError

    @abstractmethod
    def surface_area(self, dem: str, elevation: float) -> float:
        raise NotImplementedError

    @abstractmethod
    def min_max_elevations(self, dem: str) -> tuple[float, float]:
        raise NotImplementedError

    @abstractmethod
    def array(self, dem: str) -> np.ndarray:
        raise NotImplementedError

    @abstractmethod
    def histogram(self, dem: str) -> HistData:
        raise NotImplementedError

    @abstractmethod
    def create_contours(self, dem: str, polylines: str, contours: list[float]) -> None:
        raise NotImplementedError
