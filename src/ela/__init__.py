import typing

from .geoprocessing import Geoprocessor, HistData

_backend: typing.Optional[Geoprocessor] = None
"""Pick one!"""


def set_backend(be: Geoprocessor) -> None:
    global _backend
    # print("set_backend", be)
    _backend = be


def get_backend() -> typing.Optional[Geoprocessor]:
    global _backend
    # print("get_backend", _backend)
    return _backend
