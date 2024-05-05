from abc import ABC, abstractmethod
from typing import Literal

import numpy as np


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
    def histogram(self, dem: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        raise NotImplementedError
