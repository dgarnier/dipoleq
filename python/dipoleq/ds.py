from abc import ABC, abstractmethod
from dataclasses import dataclass
from imas.ids_base import IDSBase
from omas import ODS
from typing import List
import numpy as np


class DS(ABC):
    def __getitem__(self, key):
        return self._getitem(self._separate_key(key))

    def __setitem__(self, key, value):
        self._setitem(self._separate_key(key), value)

    def __getattr__(self, name):
        if name.startswith('_') and hasattr(super(), '__getattr__'):
            return super().__getattr__(name)
        return self[name]

    def __setattr__(self, name, value):
        if name.startswith('_'):
            super().__setattr__(name, value)
        else:
            self[name] = value
    
    def _separate_key(self, key):
        if isinstance(key, int):
            return [key]
        split_key = key.replace('/', '.').replace(']', '').replace('[', '.').split('.')
        return [int(k) if k.isdecimal() else k for k in split_key]
    
    @abstractmethod
    def _getitem(self, key):
        ...

    @abstractmethod
    def _setitem(self, key, value):
        ...

    @abstractmethod
    def __len__(self):
        ...
    
    @property
    @abstractmethod
    def inner(self):
        ...


@dataclass
class ImasDS(DS):
    _ids: IDSBase

    def _getitem(self, key):
        return ImasDS(self._ids)
    
    def _setitem(self, key, value):
        pass
    
    def __len__(self) -> int:
        return 0
    
    @property
    def inner(self) -> IDSBase:
        return self._ids


@dataclass
class OmasDS(DS):
    _ods: ODS

    def _getitem(self, key: List[str | int]) -> 'OmasDS':
        return self.__class__(self._ods[self._join_key(key)])
    
    def _setitem(self, key: List[str | int], value):
        assert len(key) > 0
        if len(key) == 2 and isinstance(key[1], int):
            # Deal with time-dependent arrays
            orig_value = []
            if key[0] in self._ods:
                orig_value = np.atleast_1d(self._ods[key[0]]).toList()
            if key[1] < len(orig_value):
                orig_value[key[1]] = value
            elif key[1] == len(orig_value):
                orig_value += [value]
            else:
                raise IndexError
            self._ods[key[0]] = np.atleast_1d(orig_value)
        elif len(key) == 1:
            self._ods[key[0]] = value
        else:
            OmasDS(self._ods[key[0]])._setitem(key[1:], value)
    
    def _join_key(self, key: List[str | int]) -> str:
        assert all(isinstance(k, str) or isinstance(k, int) for k in key)
        return '.'.join([str(k) for k in key])
    
    def __len__(self) -> int:
        return len(self._ods)
    
    @property
    def inner(self) -> ODS:
        return self._ods
