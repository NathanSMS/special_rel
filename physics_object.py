from four_vector import FourVector
from math import sqrt
import numpy as np

class PhysicsObject:
    c = 299_792_458  # Speed of light in meters per second

    @staticmethod
    def LorentzFactor(v, c=c) -> float:
        # Gamma = 1/sqrt(1-v^2/c^2)
        beta = v/c # Percentage of the speed of light
        radicand = 1 - beta**2 
        return sqrt(radicand, -1/2)

    def LorentzTransform():
        ...

    def __init__(self) -> None:
        pass




def main() -> None:
    ...

if __name__ == '__main__':
    main()