from abc import ABC, abstractmethod
from dataclasses import dataclass
from imas.ids_base import IDSBase
from omas import ODS


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
        split_key = key.replace('/', '.').replace(']', '').replace('[', '.').split('.')
        return [int(k) if k.isdecimal() else k for k in split_key]
    
    @abstractmethod
    def _getitem(self, key):
        ...

    @abstractmethod
    def _setitem(self, key, value):
        ...


@dataclass
class ImasDS(DS):
    _ids: IDSBase

    def _getitem(self, key):
        print(key)
        return super()._getitem(key)
    
    def _setitem(self, key, value):
        print(key, value)


@dataclass
class OmasDS(DS):
    _ods: ODS

    def _getitem(self, key) -> 'OmasDS':
        return self.__class__(self._ods[self._join_key(key)])
    
    def _setitem(self, key, value):
        self._ods[self._join_key(key)] = value
    
    def _join_key(self, key):
        assert all(isinstance(k, str) or isinstance(k, int) for k in key)
        return '.'.join([str(k) for k in key])
