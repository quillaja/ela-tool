import typing

from .geoprocessing import Geoprocessor, HistData
from .methods import (AABR, AABR2D, AAR, AAR2D, ELA, AABROriginal, AAROriginal,
                      ELAMethod)
from .slice import Slice
from .slice import create_slices_threadpool as make_slices
