from __future__ import annotations

from typing import TypeVar, Generic

O = TypeVar('O', bound='Object')
M = TypeVar('M', bound='Morphism')

''''''''''''''''''
''' CATEGORIES '''
''''''''''''''''''
'''' OBJECTS '''''
class Object:
    pass

''''''''''''''''''
'''' MORPHISMS '''
class Morphism(Generic[O]): # Object O
    def __init__(self, X : O, Y : O) -> None:
        self.X, self.Y = X, Y
    def __call__(self, x):
        pass
