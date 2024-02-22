from __future__ import annotations
from scipy.optimize import fsolve
from four_vector import FourVector
from math import sqrt
from typing import List
import numpy as np

class PhysicsObject:
    C = 299_792_458  # Speed of light in meters per second
    physics_objects = []

    def __init__(self, fourPos: FourVector = FourVector(0,0,0,0), fourVel: FourVector = FourVector(0,0,0,0), fourAcc: FourVector = FourVector(0,0,0,0), mass: float = 1, isDynamic = True):
        self.FourPos = fourPos
        self.FourVel = FourVector(0,0,0,0) if not isDynamic else fourVel
        self.FourAcc = FourVector(0,0,0,0) if not isDynamic else fourAcc
        self.mass = mass
        self.isDynamic = isDynamic
        
        self.physics_objects.append(self)

    def collide(self, other: PhysicsObject, elasticity: float = 0.9):
        if other.__class__ != self.__class__:
            raise TypeError
        else:
            if other.isDynamic:
                # Handle collision by solving 4-momentum conservation and energy conservation equations
                # P initial = P final - yields 4 simultaneous equations but 8 unknowns
                # Added constraints of |P| = |m_0*U| = m_0*c gives 2 more equations (one for each particle)
                # elasticity*KE initial = KE final ?
                pass
            else:
                # Not sure what I want to do yet, will definitely place self to the nearest point outside of other
                # Will either stop self in its tracks or reflect self's spatial velocity across the perpendicular of the line between selfs position before and after collision 
                pass

    def apply_force(self, FourForce):
        self.FourAcc = self.FourAcc + FourForce/self.mass if self.isDynamic else FourVector(0,0,0,0)

    def _time_step(self, d_ctau: float) -> None:
        if self.isDynamic:
            self.FourPos += self.FourVel*d_ctau
            self.FourVel += self.FourAcc*d_ctau
            self.FourAcc = FourVector(0,0,0,0)
        else:
            # Only add time change to position
            self.FourPos = self.FourPos + FourVector(d_ctau,0,0,0)

    def LorentzTranform(self, vectors: List[FourVector], isPosition: bool = True):
        """Transforms list of vectors to coordinates expressed in self's frame which is centered at self's position and otherwise aligned with the global frame_summary_

        Args:
            vector (FourVector): The vector to be transformed into 
            isPosition (bool, optional): _description_. Defaults to True.

        Returns:
            _type_: _description_
        """
        output = []

        # Get info for 1D lorentz transform
        c = self.FourVel.C

        gamma = self.FourVel.ct
        vx = self.FourVel.x/gamma
        vy = self.FourVel.y/gamma
        vz = self.FourVel.z/gamma
        v = sqrt(vx**2 + vy**2 + vz**2)

        row1 = [gamma,                  -gamma*vx/c,            -gamma*vy/c,            -gamma*vz/c]
        row2 = [-gamma*vx/c, 1+(gamma-1)*vx**2/v**2,    (gamma-1)*vy*vx/v*2,    (gamma-1)*vz*vx/v*2]
        row3 = [-gamma*vy/c,    (gamma-1)*vy*vx/v*2, 1+(gamma-1)*vy**2/v**2,    (gamma-1)*vz*vx/v*2]
        row4 = [-gamma*vz/c,    (gamma-1)*vz*vx/v*2,    (gamma-1)*vz*vy/v*2, 1+(gamma-1)*vz**2/v**2]
        transformation = np.asmatrix([row1,row2,row3,row4])

        for vector in vectors:
            
            if isPosition:
                # Shift vector to origin of self's coordinate system before applying transformation
                vector = vector - self.FourPos

            temp_array = np.array([vector.ct, vector.x, vector.y, vector.z], dtype=float)
            temp_array = np.reshape(temp_array, (4,1))
            transformed = np.matmul(transformation,temp_array)
            output.append(FourVector(transformed[0,0], transformed[1,0], transformed[2,0], transformed[3,0]))
        
        return output

    @classmethod
    def apply_timestep(cls, d_ctau: float):
        
        for obj in cls.physics_objects:
            obj._time_step(d_ctau)


def main() -> None:
    c = PhysicsObject.C
    test_obj = PhysicsObject(fourVel = FourVector.fourVel_from_spatialVel(0.95*c, 0, 0))
    test2_obj = PhysicsObject(fourVel = FourVector.fourVel_from_spatialVel(0.96*c, 0, 0))

    sim_time = 50
    ctau = 0.0
    d_ctau = 1
    while ctau < sim_time:
        PhysicsObject.apply_timestep(d_ctau)
        ctau += d_ctau

        test_pos_vect = test_obj.LorentzTranform(vectors=[test2_obj.FourPos])
        print(test_pos_vect[0])

        




if __name__ == '__main__':
    main()