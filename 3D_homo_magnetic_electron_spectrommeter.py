import numpy as np
import matplotlib
import scipy
import math
import importlib

IOtool = importlib.import_module("IOtool")

# in Si unit
BLX = 0.1 #magnetic field box
BLY = 0.03
BLZ = 0.02
LX = 1 #global boundary
LY = 1
LZ = 1
boxgeo = np.array([0., BLX,-0.5*BLY,0.5*BLY,-0.5*BLZ,0.5*BLZ], dtype = float) # X max, X min, Y max, Ymin, Z max, Z min
B_strength = 0.5 #Tesla
dtinbox = pow(10,-11)
dt = dtinbox
dtoutbox = pow(10,-7)
c2 = 9*pow(10,16)
m0 = 9.1*pow(10,-31)
q = 1.6*pow(10,-19)
pointsource = np.array([0.,0.,0.])
incidentEinMeV = 100.0
spectromplane = [0.0,1.0,0.,0.,0.5,0.] #a,b,c,x1,y1,y2  a(x-x1)+b(y-y1)+c(z-z1) = 0
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
    planeside = spectromplane[0]*(coordinate[0]-spectromplane[3]) + spectromplane[1]*(coordinate[1]-spectromplane[4]) + spectromplane[2]*(coordinate[2]-spectromplane[5])
    if (planeside > 0): 
        print(" out %",coordinate)
        return True
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

def updateMagneticField(coordinate,B):
    return B

def updateVelocity(state,delta_t):
    B = updateMagneticField(state.cor,B_strength)
    dpdt = np.array([-q*state.v[1]*B,-q*state.v[0]*B,0], dtype = float)
    #print(dpdt)
    p_new= momentum(state.v)+ dpdt*dt;
    v_magnitude =math.sqrt( state.v[0]*state.v[0]+state.v[1]*state.v[1]+state.v[2]*state.v[2])
    p_magnitude =math.sqrt( p_new[0]*p_new[0]+p_new[1]*p_new[1]+p_new[2]*p_new[2])
    g=p_magnitude/(v_magnitude*m0)
    #print(g)
    #print(p_new)
    state.v[0]= p_new[0]/(g*m0)
    state.v[1]= p_new[1]/(g*m0)
    state.v[2]= p_new[2]/(g*m0)
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
        print('In_box ',state.cor[0],state.cor[1],state.cor[2])
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
    print('Moveing ',state.cor[0],state.cor[1],state.cor[2])
    #print(state.cor)
    return 0 # out box but in border

def aElectronPathCalc(state,delta_t):
    pathrecord = []
    stepcounter = 0
    while(True):
        pathrecord.append(state.cor.copy())
        stepcounter += 1
        flagLocation = updateLocation(state,delta_t)
        if (flagLocation==0):
            continue
        elif (flagLocation==1):
            updateVelocity(state,delta_t)
        elif (flagLocation==2):
            print('To Spectrom and finish')
            print(' % % %',state.cor[0],state.cor[1],state.cor[2])
            break
        else:
            print('Finish')
            break
    print(stepcounter)
    return pathrecord
    
def aElectronCalc():
    state = source(pointsource, incidentEinMeV)
    pathrecord = aElectronPathCalc(state,dt)
    curve = IOtool.list_of_cor_2_3_cor_list(pathrecord)
    box = IOtool.Box()
    box.xyzmin = [0,-0.5*BLY,-0.5*BLZ]
    box.xyzmax = [BLX,0.5*BLY,0.5*BLZ]
    IOtool.plot_model(box ,curve)
    
aElectronCalc()

