from __future__ import annotations
from dataclasses import dataclass
from math import sqrt, pow

C = 299_792_458  # Speed of light in meters per second

@dataclass(frozen=True)
class FourVector:
    ct: float
    x: float
    y: float
    z: float

    @staticmethod
    def LorentzFactor(vx: float, vy: float=0, vz: float=0, c: float=C) -> float:
        """Calculates the lorentz factor for a frame moving at velocity (vx, vy, vz) with respect to the current frame

        Args:
            vx (float): x-component of velocity or full magnitude of velocity if other components are omitted
            vy (float, optional): y-component of velocity. Defaults to 0.
            vz (float, optional): y-component of velocity. Defaults to 0.
            c (float, optional): Speed of light to be used for calculation. Lower values result in more . Defaults to C.

        Returns:
            float: _description_
        """
        # Gamma = 1/sqrt(1-v^2/c^2)
        beta_squared = (vx**2 + vy**2 + vz**2)/C**2
        if beta_squared >= 1:
            raise ValueError('Magnitude of spatial velocity cannot equal or exceed the speed of light')
        radicand = 1 - beta_squared 
        return pow(radicand, -1/2)

    @classmethod
    def fourVel_from_spatialVel(cls, vx, vy, vz):
        """Returns a FourVector which is the four velocity of a frame with the spatial velocity (vx, vy, vz) as measured in the current frame

        Args:
            vx (float): x-component of measured velocity
            vy (float): y-component of measured velocity
            vz (float): z-component of measured velocity

        Returns:
            FourVector: The four velocity, U = (dct, ux, uy, uz) for a given spatial velocity (vx, vy, vz)
        """
        # dU/dTau = dt/dTau * dU/dt = gamma*[ct, vx, vy, vz]
        gamma = cls.LorentzFactor(vx, vy, vz)
        u_ct = gamma*C
        u_x = gamma*vx
        u_y = gamma*vy
        u_z = gamma*vz
        return FourVector(u_ct, u_x, u_y, u_z)

    def square_mag(self) -> float:
        """Calculates four-vector magnitude squared using Minkowski Metric with mostly minuses convention

        Returns:
            float: Square magnitude of vector under minkowski mostly minuses convention
        """
        return self.ct**2 - self.x**2 - self.y**2 - self.z**2
    
    def spatialMag(self) -> float:
        """Returns the spatial magnitude of a four vector

        Returns:
            float: Spatial magnitude of vector
        """
        return sqrt(self.x**2 + self.y**2 + self.z**2)
    
    def __abs__(self) -> float:
        """Returns absolute value of magnitude of vector using Minkowski Mostly Minuses convention, |ct^2 - x^2 - y^2 - z^2|.
        If self is a FourVelocity vector, this should always be C

        Returns:
            float: Absolute vector magnitude
        """
        return sqrt(abs(self.square_mag()))
    
    def __rmul__(self, other: float|int) -> FourVector:
        """Returns four vector which is self scaled by other

        Args:
            other (float | int): Scalar to multiply four vector by

        Returns:
            FourVector: result of scalar multiplication other*self
        """
        if other.__class__ in (float, int):
            return FourVector(self.ct*other, self.x*other, self.y*other, self.z*other)
        else:
            return NotImplemented
    
    def dot(self, other: FourVector) -> float:
        """Calculates dor product between self and other using Minkowski metric and mostly minuses convention

        Args:
            other (FourVector): FourVector to be dotted with self

        Returns:
            float: Dot produt of self and other using minkowski metric and mostly minuses convention
        """
        if self.__class__ == other.__class__:
            return self.ct*other.ct - self.x*other.x - self.y*other.y - self.z*other.z
        else:
            return NotImplemented

    def __add__(self, other: FourVector) -> FourVector:
        """Returns sum of 2 FourVectors

        Args:
            other (FourVector): FourVector to be added to self

        Returns:
            FourVector: element-wise sum of self and other
        """
        if self.__class__ == other.__class__:
            return FourVector(self.ct+other.ct, self.x+other.x, self.y+other.y, self.z+other.z)
        else:
            return NotImplemented
        
    def __sub__(self, other) -> FourVector:
        """Returns difference of 2 FourVectors

        Args:
            other (FourVector): FourVector to be subtracted from self

        Returns:
            FourVector: element-wise subtraction of other from self
        """
        return self + -1*other

def main() -> None:
    test1 = FourVector.fourVel_from_spatialVel(vx = 0.1*C, vy=0.1*C, vz=0.1*C)
    print(f'Test: {test1}')
    print(f'Magnitude: {abs(test1)}')
    print(f'Speed of light: {C}')

if __name__ == '__main__':
    main()