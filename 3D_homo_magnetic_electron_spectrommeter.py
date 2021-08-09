import numpy as np
import matplotlib
import scipy
import math

# in Si unit
BLX = 0.1
BLY = 0.05
BLZ = 0.05
LX = 1
LY = 1
LZ = 1
boxgeo = np.array([0., BLX,-0.5*BLY,0.5*BLY,-0.5*BLZ,0.5*BLZ], dtype = float) # X max, X min, Y max, Ymin, Z max, Z min
B = 1 #Tesla
dtinbox = pow(10,-11)
dt = dtinbox
dtoutbox = pow(10,-7)
c2 = 9*pow(10,16)
m0 = 9.1*pow(10,-31)
q = 1.6*pow(10,-19)
# Location 'noFieldRegion', 'FieldRegion', 'Spectrometer', 'Outboder', 'EnergyTooLow'

class particle_state():
    def __init__(self):
        self.cor=np.array([0.,0.,0.], dtype = float)
        self.v=np.array([0.,0.,0.], dtype = float)

def IsInMagBox(coordinate):
    if (coordinate[0]<boxgeo[0] or coordinate[0]>boxgeo[1]): return False
    if (coordinate[1]<boxgeo[2] or coordinate[1]>boxgeo[3]): return False
    if (coordinate[2]<boxgeo[4] or coordinate[2]>boxgeo[5]): return False
    return True

def ToSpectrom(coordinate):
    if (coordinate[1] < -0.05): return True
    return False

def IsOutBoder(coordinate):
    if (coordinate[0]<-0.5*LX or coordinate[0]>0.5*LX): return True
    if (coordinate[1]<-0.5*LY or coordinate[1]>0.5*LY): return True
    if (coordinate[2]<-0.5*LZ or coordinate[2]>0.5*LZ): return True
    return False

def gamma(velocity):
    magnitude_square = velocity[0]*velocity[0]+velocity[1]*velocity[1]+velocity[2]*velocity[2]
    return 1./math.sqrt(1.-(magnitude_square/c2))

def momentum(velocity):
    g = gamma(velocity)
    gm0 = g*m0
    return np.array([gm0*velocity[0],gm0*velocity[1],gm0*velocity[2]], dtype = float)

def updateVelocity(velocity,delta_t):
    dpdt = np.array([-q*velocity[1]*B,-q*velocity[0]*B,0], dtype = float)
    #print(dpdt)
    p_new= momentum(velocity)+ dpdt*dt;
    v_magnitude =math.sqrt( velocity[0]*velocity[0]+velocity[1]*velocity[1]+velocity[2]*velocity[2])
    p_magnitude =math.sqrt( p_new[0]*p_new[0]+p_new[1]*p_new[1]+p_new[2]*p_new[2])
    g=p_magnitude/(v_magnitude*m0)
    #print(g)
    #print(p_new)
    velocity[0]= p_new[0]/(g*m0)
    velocity[1]= p_new[1]/(g*m0)
    velocity[2]= p_new[2]/(g*m0)
    #print(velocity)
    
def calcElectronVel(energy_in_MeV):
    mc2 = 0.511 #Mev/c^2
    return math.sqrt(c2-c2/math.pow((1+energy_in_MeV/mc2),2))
    
def source(coordinate, energy_in_MeV, thetaY = 0., thetaZ = 0.):
    state = particle_state()
    state.cor = coordinate
    state.v[0] =  calcElectronVel(energy_in_MeV)*math.cos(thetaY)*math.cos(thetaZ)
    state.v[1] =  calcElectronVel(energy_in_MeV)*math.sin(thetaY)*math.cos(thetaZ)
    state.v[2] =  calcElectronVel(energy_in_MeV)*math.cos(thetaY)*math.sin(thetaZ)
    return state

def updateLocation(state,delta_t):
    state.cor[0]+=state.v[0]*delta_t
    state.cor[1]+=state.v[1]*delta_t
    state.cor[2]+=state.v[2]*delta_t
    if (IsInMagBox(state.cor)): 
        print('In box '+str(state.cor))
        return 1
    if (ToSpectrom(state.cor)): 
        print('To Spectrom')
        return 2
    if (IsOutBoder(state.cor)): 
        print('Out of boder') 
        return 3
    if (state.v[0]<0.01): 
        print('V too low') 
        return 4
    print(state.cor)
    return 0 # out box but in border

def aElectronPathCalc(state,delta_t):
    stepcounter = 0
    while(True):
        stepcounter += 1
        flagLocation = updateLocation(state,delta_t)
        if (flagLocation==0):
            dt = dtoutbox
            continue
        elif (flagLocation==1):
            updateVelocity(state.v,delta_t)
        elif (flagLocation==2):
            print('To Spectrom and finish')
            print(state.cor)
            break
        else:
            print('Finish')
            break
    print(stepcounter)
    
def aElectronCalc():
    state = source(np.array([0.,0.,0.], dtype = float),100)
    aElectronPathCalc(state,dt)
    
aElectronCalc()

