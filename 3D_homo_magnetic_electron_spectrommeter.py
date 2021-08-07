import numpy as np
import matplotlob
import scipy
import math

# in Si unit
boxgeo = np.array([0., 0.2,0.,0.1,0.,0.2], dtype = float) # X max, X min, Y max, Ymin, Z max, Z min
B = 10. #in Z direction
dt = pow(10,-9)

def particle_state():
    def __init__(self):
        self.cor=np.array([0.,0.,0.], dtype = float)
        self.v=np.array([0.,0.,0.], dtype = float)

def IsInMagBox(coordinate):
    if (coordinate[0]<boxgeo[0] or coordinate[0]>boxgeo[1]): return False
    if (coordinate[1]<boxgeo[2] or coordinate[1]>boxgeo[3]): return False
    if (coordinate[2]<boxgeo[4] or coordinate[2]>boxgeo[5]): return False
    return True

def ToSpectrom(coordinate):
    if (coordinate[0]<boxgeo[0] or coordinate[0]>boxgeo[1]): return False
    if (coordinate[1]<boxgeo[2] or coordinate[1]>boxgeo[3]): return False
    if (coordinate[2]<boxgeo[4] or coordinate[2]>boxgeo[5]): return False
    return True


def gamma(velocity):
    magnitude_square = velocity[0]*velocity[0]+velocity[1]*velocity[1]+velocity[2]*velocity[2]
    c_square = 9*pow(10,16)
    return 1/math.sqrt(1-(magnitude_square/c_square))

def UpdateVelocity(velocity,delta_t):
    q = 1.6*pow(10,-19)
    force = q /pow(10,16)*velocity[0]*B #Lorentz force
    g = gamma(velocity)
    m = 9.1*pow(10,-31)
    acceleration = (1/(m*g))*(force/g)
    velocity[1]=velocity[1]+0.5*delta_t*acceleration #middle velocity
    
def source(coordinate, energy_in_MeV, thetaX = 0., thetaY = 0.):
    state = particle_state()
    return state

def updatelocation(state):
    if (IsInMagBox(state.cor)): 
        print('In box')
        UpdateVelocity(state.v,dt)
    if (ToSpectrom(state.cor)): print('To Spectrom')
    if (state.v<0.01): 
        print('V too low') 
        return -1
    print('update cor')
    return 0

