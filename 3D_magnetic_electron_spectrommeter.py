import numpy as np
import matplotlib
import scipy
import math
import importlib

Plottool = importlib.import_module("Plot_tool_magnetic_spectrommeter")
MagneticField = importlib.import_module("magnetic_field")

# in Si unit
BLX = 0.01 #magnetic field box
BLY = 0.04
BLZ = 0.02
LX = 1 #global boundary
LY = 1
LZ = 1
boxgeo = np.array([0., BLX,-0.5*BLY,0.5*BLY,-0.5*BLZ,0.5*BLZ], dtype = float) # X max, X min, Y max, Ymin, Z max, Z min
spectromplane = [1,0,0.,0.17,0.0,0.] #a,b,c,x1,y1,y2  a(x-x1)+b(y-y1)+c(z-z1) = 0

dtinbox = pow(10,-14)
dt = dtinbox
dtoutbox = pow(10,-13)


pointsource = np.array([0.,0.,0.])
incidentEinMeV = 100.0

B_strength = 0.22 #Tesla
Is_magnetic_homo = True
MagneticField2D = []
MagneticFilePath = 'magnet_1_measurement_for_test.csv'
gp0 = 0.;
ggamma0 = 0.;
gv0 = 0;

Is2D = True
# Location 'noFieldRegion', 'FieldRegion', 'Spectrometer', 'Outboder', 'EnergyTooLow'

# Physical constant
c2 = 2.9979*2.9979*pow(10,16)
m0 = 9.1*pow(10,-31)
q = 1.6*pow(10,-19)

class particle_state():
    def __init__(self,cornum = None, vnum = None):
        if vnum is None:
            self.cor=np.array([0.,0.,0.], dtype = float)
            self.v=np.array([0.,0.,0.], dtype = float)
        else:
            self.cor=cornum
            self.v=vnum

def IsInMagBox(coordinate):
    if (coordinate[0]<boxgeo[0] or coordinate[0]>boxgeo[1]): return False
    if (coordinate[1]<boxgeo[2] or coordinate[1]>boxgeo[3]): return False
    if (coordinate[2]<boxgeo[4] or coordinate[2]>boxgeo[5]): return False
    return True

def ToSpectrom(coordinate):
    planeside = spectromplane[0]*(coordinate[0]-spectromplane[3]) + spectromplane[1]*(coordinate[1]-spectromplane[4]) + spectromplane[2]*(coordinate[2]-spectromplane[5])
    if (planeside > 0): 
        print(" out ",coordinate)
        return True
    return False

def IsOutBoder(coordinate):
    if (coordinate[0]<-0.5*LX or coordinate[0]>0.5*LX): return True
    if (coordinate[1]<-0.5*LY or coordinate[1]>0.5*LY): return True
    if (coordinate[2]<-0.5*LZ or coordinate[2]>0.5*LZ): return True
    return False

def gamma_vec(velocity):
    magnitude_square = velocity[0]*velocity[0]+velocity[1]*velocity[1]+velocity[2]*velocity[2]
    return 1./math.sqrt(1.-(magnitude_square/c2))

def gamma_pure(velocity):
    return 1./math.sqrt(1.-(velocity*velocity/c2))

def momentum(velocity):
    g = gamma_vec(velocity)
    gm0 = g*m0
    return np.array([gm0*velocity[0],gm0*velocity[1],gm0*velocity[2]], dtype = float)

def updateMagneticField(coordinate,B):
    if (Is_magnetic_homo and len(MagneticField2D)>0):
        return B
    fieldX_ele = len(MagneticField2D[0])
    fieldY_ele = len(MagneticField2D)
    gridx = (int)((coordinate[0]/BLX)*fieldX_ele)
    gridy = (int)((coordinate[1]/BLY)*fieldY_ele)
    #print("B field to be % % %",gridx, gridy, float(MagneticField2D[gridx][gridy]) )
    return float(MagneticField2D[gridx][gridy])
    

def updateVelocity(state,delta_t):
    B = updateMagneticField(state.cor,B_strength)
    dpdt = np.array([-q*state.v[1]*B,-q*state.v[0]*B,0], dtype = float)
    #print(dpdt)
    p_new= momentum(state.v)+ dpdt*dt;
    p_magnitude =math.sqrt( p_new[0]*p_new[0]+p_new[1]*p_new[1]+p_new[2]*p_new[2])
    v_magnitude =math.sqrt( state.v[0]*state.v[0]+state.v[1]*state.v[1]+state.v[2]*state.v[2])
    g=p_magnitude/(v_magnitude*m0)  
    state.v[0]= p_new[0]/(g*m0)
    state.v[1]= p_new[1]/(g*m0)
    state.v[2]= p_new[2]/(g*m0)
    
def updateVelocity2(state,delta_t):
    B = updateMagneticField(state.cor,B_strength)
    g = gamma_vec(state.v)
    state.v[0] = state.v[0] + q/g/m0*state.v[1]*delta_t*B
    state.v[1] = state.v[1] - q/g/m0*state.v[0]*delta_t*B  

    
def calcElectronVel(energy_in_MeV):
    mc2 = 0.511 #Mev/c^2
    return math.sqrt(1-1/math.pow((1+energy_in_MeV/mc2),2))*2.9979*pow(10,8);
    
def source(coordinate, energy_in_MeV, thetaY = 0., thetaZ = 0.):
    state = particle_state()
    state.cor = coordinate
    state.v[0] =  calcElectronVel(energy_in_MeV)*math.cos(thetaY)*math.cos(thetaZ)
    state.v[1] =  calcElectronVel(energy_in_MeV)*math.sin(thetaY)*math.cos(thetaZ)
    state.v[2] =  calcElectronVel(energy_in_MeV)*math.cos(thetaY)*math.sin(thetaZ)
    global ggamma0 
    ggamma0 = gamma_vec(state.v)
    return state

def check_pos(state: particle_state):
    if (IsInMagBox(state.cor)): 
        #print('In_box ',state.cor[0],state.cor[1],state.cor[2])
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
    return 0 # out box but in border

def updateLocation(state,delta_t):
    state.cor[0]+=state.v[0]*delta_t
    state.cor[1]+=state.v[1]*delta_t
    state.cor[2]+=state.v[2]*delta_t
    return check_pos(state)

def updateLocation4thRK(state,delta_t):
    global ggamma0
    k1x = delta_t * (state.v)
    k1p = delta_t * (q * np.cross(state.v, np.array([0.,0.,updateMagneticField(state.cor,B_strength)], dtype = float)))
    vtmp = state.v + k1p/2/ggamma0/m0
    xtmp = state.cor + k1x/2
    k2x = delta_t * (vtmp)
    k2p = delta_t * (q * np.cross(vtmp, np.array([0.,0.,updateMagneticField(xtmp,B_strength)], dtype = float)))
    vtmp = state.v + k2p/2/ggamma0/m0
    xtmp = state.cor + k2x/2
    k3x = delta_t * (vtmp)
    k3p = delta_t * (q * np.cross(vtmp, np.array([0.,0.,updateMagneticField(xtmp,B_strength)], dtype = float)))
    vtmp = state.v + k3p/ggamma0/m0
    xtmp = state.cor + k3x
    k4x = delta_t * (vtmp)
    k4p = delta_t * (q * np.cross(vtmp, np.array([0.,0.,updateMagneticField(xtmp,B_strength)], dtype = float)))
    state.cor = state.cor + k1x/6 + k2x/3 + k3x/3 + k4x/6
    state.v = state.v + (k1p/6 + k2p/3 + k3p/3 + k4p/6)/ggamma0/m0
    return check_pos(state)

def aElectronPathCalc(state,delta_t):
    pathrecord = []
    stepcounter = 0
    while(True):
        pathrecord.append(state.cor.copy())
        stepcounter += 1
        flagLocation = check_pos(state);
        if (flagLocation==0):
            updateLocation(state,delta_t)
            continue
        elif (flagLocation==1):
            updateLocation4thRK(state,delta_t)
        elif (flagLocation==2):
            print('To Spectrom and finish')
            print(state.cor[0],state.cor[1],state.cor[2])
            break
        else:
            print(state.cor[0],state.cor[1],state.cor[2])
            print('Finish')
            break
    print(stepcounter)
    return pathrecord
    
def aElectronCalc():
    electronpaths = []
    count = 1
    energies = []
    for i in range (1,11,1):
        state = source(pointsource.copy(), i)
        #print(state.v[0])
        pathrecord = aElectronPathCalc(state,dt)
        curve = Plottool.list_of_cor_2_3_cor_list(pathrecord)
        electronpaths.append(curve)
        energies.append(count)
        count+=1
    box = Plottool.Box()
    box.xyzmin = [0,-0.5*BLY,-0.5*BLZ]
    box.xyzmax = [BLX,0.5*BLY,0.5*BLZ]
    if (Is2D):
        Plottool.plot_model2D(box ,electronpaths)
        Plottool.plot_energy_spectrom(electronpaths, energies)
    else:
        Plottool.plot_model3D(box ,electronpaths)
        
        
MagneticField2D = MagneticField.readMagnitCsv(MagneticFilePath)
aElectronCalc()

