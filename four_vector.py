from __future__ import annotations
from dataclasses import dataclass
from math import sqrt

@dataclass(frozen=True)
class FourVector:
    ct: float
    x: float
    y: float
    z: float

    def square_mag(self) -> float:
        return self.ct**2 - self.x**2 - self.y**2 - self.z**2
    
    def __abs__(self) -> float:
        return sqrt(abs(self.square_mag()))
    
    def __rmul__(self, other) -> FourVector:
        if other.__class__ in (float, int):
            return FourVector(self.ct*other, self.x*other, self.y*other, self.z*other)
        else:
            return NotImplemented
    
    def dot(self, other: FourVector) -> float:
        if self.__class__ == other.__class__:
            return self.ct*other.ct - self.x*other.x - self.y*other.y - self.z*other.z
        else:
            return NotImplemented

    def __add__(self, other) -> FourVector:
        if self.__class__ == other.__class__:
            return FourVector(self.ct+other.ct, self.x+other.x, self.y+other.y, self.z+other.z)
        else:
            return NotImplemented
        
    def __sub__(self, other) -> FourVector:
        return self + -1*other



def main() -> None:
    test1 = FourVector(1,2,3,4)
    test2 = FourVector(1,2,3,4)
    print(test2.dot(test1))

if __name__ == '__main__':
    main()