from four_vector import FourVector
from math import sqrt
import numpy as np

class PhysicsObject:
    C = 299_792_458  # Speed of light in meters per second
    dynamic_objects = []
    static_objects = []

    def __init__(self, fourPos: FourVector = FourVector(0,0,0,0), fourVel: FourVector = FourVector(0,0,0,0), fourAcc: FourVector = FourVector(0,0,0,0), mass: float = 1, isDynamic = True):
        self.FourPos = fourPos
        self.FourVel = fourVel
        self.FourAcc = fourAcc
        self.mass = mass
        self.isDynamic = isDynamic
        if isDynamic:
            self.dynamic_objects.append(self)
        else:
            self.static_objects.append(self)

    def apply_force(self, FourForce):
        self.FourAcc = self.FourAcc + FourForce/self.mass

    def _time_step(self, d_ctau: float) -> None:
        if self.isDynamic:
            self.FourPos += self.FourVel*d_ctau
            self.FourVel += self.FourAcc*d_ctau
            self.FourAcc = 0
            self.animationCount += 1

    def LorentzTranform(self, vector, c=C):
        # Get info for 1D lorentz transform
        v = self.FourVel.spatialMag()
        gamma = self.FourVel.LorentzFactor(v,c)
        beta = v/c

        vx = self.FourVel.x
        vy = self.FourVel.y
        vz = self.FourVel.z

        temp_array = np.array([vector.ct, vector.x, vector.y, vector.z], dtype=float, shape=(4, 1))

        row1 = [gamma,                  -gamma*vx/c,            -gamma*vy/c,            -gamma*vz/c]
        row2 = [-gamma*vx/c, 1+(gamma-1)*vx**2/v**2,    (gamma-1)*vy*vx/v*2,    (gamma-1)*vz*vx/v*2]
        row3 = [-gamma*vy/c,    (gamma-1)*vy*vx/v*2, 1+(gamma-1)*vy**2/v**2,    (gamma-1)*vz*vx/v*2]
        row4 = [-gamma*vz/c,    (gamma-1)*vz*vx/v*2,    (gamma-1)*vz*vy/v*2, 1+(gamma-1)*vz**2/v**2]
        transformation = np.asmatrix([row1,row2,row3,row4])
        transformed = np.matmul(transformation,temp_array)
        return FourVector(transformed[0], transformed[1], transformed[2], transformed[3])

    @classmethod
    def apply_timestep(cls, d_ctau: float):
        
        for obj in cls.dynamic_objects:
            # obj._time_step()
            ...



def main() -> None:
    c = PhysicsObject.C
    test_obj = PhysicsObject(fourVel = FourVector.fourVel_from_spatialVel(0.1*c, 0.1*c, 0.1*c))
    print(test_obj.FourVel)
    print(test_obj.FourVel.square_mag())

if __name__ == '__main__':
    main()