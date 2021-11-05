#Detailed Disturbance Torque Calculations
import Attitude_Kinematics as att
import numpy as np
import datetime
import csv

class Orientation:
    def __init__(self,  date            = datetime.datetime(1,1,1,1,1,1),
                        position        = np.matrix([ [0], [0], [0] ]), 
                        velocity        = np.matrix([ [0], [0], [0] ]), 
                        euler_angles    = np.array([0,0,0])):
        self.date           = date
        self.position       = position
        self.velocity       = velocity
        self.euler_angles   = euler_angles

class Plate:
    def __init__(self,  area            = 0.0,
                        normal          = np.matrix([ [1], [0], [0] ]), 
                        lever           = np.matrix([ [1], [0], [0] ]) , 
                        r_coef          = np.array([0,0,0])):
        self.area   = area
        self.normal = normal
        self.lever  = lever
        self.r_coef = r_coef


#importing atmospheric data from csv
atmos_table = []
with open('Atmospheric_density_model.csv', 'r', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in reader:
        atmos_table.append([ float(row[0]), float(row[1]), float(row[2]) ])


#constants
PHI = 1361 #solar constant TODO units
c = 3*10**8 #speed of light (m/s)
Wo = np.matrix([ [0], [0], [0.000072921158553]]) #earths rotation in rads/s
Cd = 2.25 #typical coefficient of drag ~2-2.5

def sind(deg):
    return np.sin(np.deg2rad(deg))


def cosd(deg):
    return np.cos(np.deg2rad(deg))


#gives sattelite-sun vector in ECI frame
def get_rsato(date,r_sattoearth):
    Tut1    = (att.Get_Julian_Datetime(date) - 2451545) / 36525
    Mo      = 357.5277233 + 35999.05034*Tut1  #all of the subsequent calculations are in terms of degrees
    Ep      = 23.439291 - 0.0130042*Tut1
    Phi_o   = 280.46 + 36000.771*Tut1
    Phi_ecl = Phi_o + 1.914666471*sind(Mo) + 0.019994643*sind(2*Mo)

    roo_unitvec = np.matrix([ [cosd(Phi_ecl)], [cosd(Ep)*sind(Phi_ecl)], [sind(Ep)*sind(Phi_ecl)] ])
    roo_mag = 1.000140612 - 0.016708617*cosd(Mo) - 0.000139589*cosd(2*Mo)

    r_sattoearth = np.multiply(1/(1.496*10**8), r_sattoearth) #r_sattoearth given in km converted to AU
    r_sato = np.multiply(roo_mag, roo_unitvec) - r_sattoearth

    return r_sato


#angle_vec is angles arranged in the order which they are applied, order_vec is the order of axes for the transformation
def get_dcm(angle_vec, order_vec): 
    #this functions uses radians... degrees must be converted before input 
    ang = np.array([0.0,0.0,0.0]) #if these are floats then float inputs dont work?
    for i in range(3):
        ang[i] = angle_vec[order_vec[i]-1]

    z_DCM = np.matrix([
        [ np.cos(ang[0]), np.sin(ang[0]), 0 ],
        [-np.sin(ang[0]), np.cos(ang[0]), 0 ],
        [ 0             , 0             , 1 ]
    ])
    y_DCM = np.matrix([
        [ np.cos(ang[1]), 0,-np.sin(ang[1]) ],
        [ 0             , 1, 0              ],
        [ np.sin(ang[1]), 0, np.cos(ang[1]) ]
    ])
    x_DCM = np.matrix([
        [  1, 0             , 0             ],
        [  0, np.cos(ang[2]), np.sin(ang[2])],
        [  0,-np.sin(ang[2]), np.cos(ang[2])]
    ])
    
    DCM_array = np.array([ x_DCM, y_DCM, z_DCM ]) # x==1, y==2, z==3 for purpose of DCM description
    DCM = np.matrix([
        [ 1, 0, 0],
        [ 0, 1, 0],
        [ 0, 0, 1]
    ])
    for i in range(3):
        DCM = np.matmul( DCM_array[order_vec[i]-1], DCM)

    return DCM


#DCM which converts orbital to inertial frame (often you want the reverse transformation so transpose)
def get_Aio(position, velocity): 
    Aio = np.matrix([
        [0,0,0],
        [0,0,0], 
        [0,0,0]
    ])
    Aio[:,0] = np.multiply( -position, 1/np.linalg.norm(position)) 
    h = np.cross(position.transpose(), velocity.transpose())
    Aio[:,1] = np.multiply( -h, 1/np.linalg.norm(h)).transpose()
    Aio[:,2] = np.cross( Aio[:,1].transpose(), Aio[:,0].transpose() ).transpose()

    return Aio


def get_density(alt):

    i=0
    while alt >= atmos_table[i,1]:
        i+=1
    
    h0 = atmos_table[i,0]
    H = atmos_table[i,2]
    rho_0 = atmos_table[i,1]

    rho = rho_0*np.exp(-(alt-h0)/H)
    return rho

#calculation of solar torque on plate in body frame
def get_solartorque(sat_orientation, plate):
    #transformation Abi
    angles = sat_orientation.euler_angles
    angles = np.multiply(np.pi/180, angles)
    order = np.array([3,2,1])
    Abo = get_dcm(angles, order)
    Aoi = get_Aio(sat_orientation.position, sat_orientation.velocity).transpose()
    Abi = np.matmul( Abo, Aoi)
    
    #sun-sat vector in body frame
    r_sato = get_rsato(sat_orientation.date, sat_orientation.position)
    s = np.matmul( Abi, r_sato)

    #incidence angle of sun-sat vector and plate normal
    cos_theta_srp = np.matmul( plate.normal.transpose(), s )[0,0]

    #calculation of solar pressure force of plate (in body frame)
    r_diff = plate.r_coef[0]#note order
    r_spec = plate.r_coef[1]
    F_srp = -PHI/(c*np.linalg.norm(r_sato))*plate.area*( 2*(r_diff/3+r_spec*cos_theta_srp) + (1-r_spec)*s )*max(cos_theta_srp,0)

    #calculation of solar presure torque of plate
    L_srp = np.cross( plate.lever.transpose(), F_srp.transpose() ).transpose()

    return L_srp


def get_aerotorque(sat_orientation, plate):
    #getting body orientation w.r.t. ECI
    Aoi = get_Aio(sat_orientation.position, sat_orientation.velocity).transpose()
    Abo = get_dcm(sat_orientation.euler_angles, np.array([3,2,1]))
    Abi = np.matmul(Abo, Aoi)
    
    #calculate Vrel in body
    v_reli = sat_orientation.position + np.cross( Wo.transpose(), sat_orientation.position.transpose() )
    v_relb = np.matmul( Abi, v_reli ) 

    #calculation force on plate due to aerodynamic pressure F_aero
    cos_theta_aero = 1/np.linalg.norm(v_relb)*np.matmul( plate.normal.tranpose(), v_relb)
    alt = np.linalg.norm(sat_orientation.position) - 6371.0
    rho = get_density(alt)
    F_aero = -0.5*rho*Cd*np.linalg.norm(v_relb)*plate.area*max(cos_theta_aero,0)*v_relb

    #calculation of torque on plate due to aerodynamic pressure L_aero
    L_aero = np.cross( sat_orientation.lever.transpose(), F_aero.transpose() ).transpose()
    return L_aero



if __name__=="__main__":

    date = datetime.datetime(2021,11,3,10,9,0)
    r_sat = np.matrix([ [6971], [0], [0]])

    r_sato = get_rsato(date,r_sat)
    r_sato = np.multiply( 1/np.linalg.norm(r_sato), r_sato)

    angles = np.array([45,45,45])
    angles = np.multiply(np.pi/180, angles)
    order = np.array([3,2,1])
    Abo = get_dcm(angles, order)

    position = np.matrix([
        [6971],
        [0],
        [0]
    ])
    velocity = np.matrix([
        [0],
        [7.5615],
        [0],
    ])
    Aoi = get_Aio(position, velocity).transpose()
    Abi = np.matmul( Abo, Aoi)
    print(Abi)

    r_sato_b = np.matmul(Abi, r_sato)
    print(r_sato_b)
    
    sat_orientation = Orientation() #this will get initialized with default values... must change to ones you want
    sat_orientation.date = date
    sat_orientation.position = position #describes position relative to center of earth in ECI frame
    sat_orientation.velocity = velocity  #these are in km and km/s units 
    sat_orientation.euler_angles = angles #in degrees
    print(sat_orientation.date)
    print(sat_orientation.euler_angles)

    plate = Plate()
    plate.area = 1
    plate.lever = 1/np.sqrt(3)*np.matrix([ [1],[1],[1] ])
    plate.normal = 1/np.sqrt(3)*np.matrix([ [1],[1],[1] ])
    plate.r_coef = np.array([ 1/3, 1/3, 1/3 ])

    print(get_solartorque(sat_orientation,plate)) #I think this works