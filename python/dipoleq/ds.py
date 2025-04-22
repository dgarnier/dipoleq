from abc import ABC, abstractmethod
from dataclasses import dataclass
from imas.ids_base import IDSBase
from imas.ids_primitive import IDSPrimitive
from omas import ODS
from typing import Any, List
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
    
    def _wrap_object(self, o: Any, force: bool) -> Any:
        if isinstance(o, ODS):
            return OmasDS(o)
        if isinstance(o, IDSPrimitive) and not force:
            return o.value
        if isinstance(o, IDSBase):
            return ImasDS(o)
        if force:
            raise TypeError(f"Unsupported type {type(o)} of {o} in force wrap mode")
        return o
    
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

    def _getitem(self, key: List[str | int], force: bool = False) -> Any:
        if len(key) == 1:
            if isinstance(key[0], int) and key[0] == len(self._ids):
                self._ids.resize(len(self._ids)+1)
            return self._wrap_object(self._ids[key[0]], force)
        return self._getitem([key[0]])._getitem(key[1:])

    def _setitem(self, key: List[str | int], value: Any) -> None:
        if len(key) == 1:
            if isinstance(key[0], int) and key[0] == len(self._ids):
                self._ids.resize(len(self._ids)+1)
            self._ids[key[0]] = value
        else:
            self._getitem([key[0]], True)._setitem(key[1:], value)

    def _join_key(self, key: List[str | int]) -> str:
        result = ''
        for i, k in enumerate(key):
            if isinstance(k, int):
                result += f'[{k}]'
            elif isinstance(k, str):
                if i != 0:
                    result += '/'
                result += k
            else:
                raise ValueError('Invalid key type')
        return result

    def __len__(self) -> int:
        return len(self._ids)

    @property
    def inner(self) -> IDSBase:
        return self._ids


@dataclass
class OmasDS(DS):
    _ods: ODS

    def _getitem(self, key: List[str | int], force: bool = False) -> 'OmasDS':
        return self._wrap_object(self._ods[self._join_key(key)], force)
    
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
            self._getitem([key[0]], True)._setitem(key[1:], value)
    
    def _join_key(self, key: List[str | int]) -> str:
        assert all(isinstance(k, str) or isinstance(k, int) for k in key)
        return '.'.join([str(k) for k in key])
    
    def __len__(self) -> int:
        return len(self._ods)
    
    @property
    def inner(self) -> ODS:
        return self._ods
