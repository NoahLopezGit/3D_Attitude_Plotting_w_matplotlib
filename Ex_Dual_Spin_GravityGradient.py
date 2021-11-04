#Numerically Solving the EOMs of a Dual Spin Satellite with a Gravity Gradient TODO: this shit broken yo
import Attitude_Kinematics as att
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



def Get_Wdot(euler_angles, w_vec, Jmat, Jmat_inv, wheel_h, w_oi, u, R):  #TODO: make a data structure for this

    #getting DCM from orbital to body frame
    A_bo = Get_DCM(euler_angles)

    #first term in exnnp.matrix(np.cross( w_bi.transpose(), H_bi.transpose())).transpose()
    Jvec = np.array([ Jmat[0,0], Jmat[1,1], Jmat[2,2] ])
    T_gravity = Get_GravityTorque(euler_angles, Jvec, u, R)

    #second term in exn
    w_bi = w_vec + np.matmul(A_bo, w_oi).transpose()
    H_bi = wheel_h + np.matmul( Jmat, w_bi)

    #3rd term in exnt_vec, euler_array[1,:], 'r',
    A_bo_dot = np.matmul( Get_Xss_Matrix(w_vec), A_bo )
    A_dot_w_oi = np.matmul( A_bo_dot, w_oi)
    hdot = T_gravity - np.matrix(np.cross( w_bi.transpose(), H_bi.transpose())).transpose() - A_dot_w_oi.transpose()
    wdot = np.matmul( Jmat_inv, hdot)
    return wdot


def Get_DCM(euler_angles):
    e = euler_angles  #dont mind this lol
    #with small angle approximation for 3-2-1 euler angles
    DCM = np.matrix([
        [ 1     , e[2,0],-e[1,0] ],
        [-e[2,0], 1     , e[0,0] ],
        [ e[1,0],-e[0,0], 1      ]
    ])
    return DCM


def Get_Xss_Matrix(vec):
    cross_product_matrix = np.matrix([
        [ 0       ,-vec[2,0] , vec[1,0] ],
        [ vec[2,0], 0        ,-vec[0,0] ],
        [-vec[1,0], vec[0,0] , 0        ]
    ])
    return cross_product_matrix


def Get_Euler_Rates(euler_vec, w_vec):
    w2rates = np.matrix([
        [ 0, euler_vec[0,0], 1             ],
        [ 0, 1             ,-euler_vec[0,0]],
        [ 1, 0             , euler_vec[1,0]]
    ])
    euler_rates = np.matmul( w2rates, w_vec )
    return euler_rates.transpose()


def Get_W_oi():
    #TODO??
    return


def Get_GravityTorque(euler_vec, J_vec, u, R):
    torque_vec = np.matrix([
        [ (J_vec[2]-J_vec[1])*euler_vec[0,0] ],
        [ (J_vec[2]-J_vec[1])*euler_vec[1,0] ],
        [ 0]
    ])
    gravity_torque = np.multiply( (3*u/R**3), torque_vec)
    return gravity_torque



if __name__=="__main__":

    u = 3.986*10**14  #earth gravitational constant
    R = (6371 + 650)*10**3  #radius of orbit
    n = np.sqrt(u/R**3) #angular velocity of orbit
    #get initial orbital velocity LVLH vs. Inertial
    w_oi = np.array([ 0, -n, 0 ])#Get_W_oi(orbital_quaternion)TODO: why are they asking for function here
    #Implement formulas to calculate w_dot body/orbital frame
    total_time = (2*np.pi/n) #seconds for 10 orbits
    timesteps = 1001 #TODO: Need 10 orbits @ 2pi/n time per orbit
    t_vec = np.matrix(np.linspace(0,total_time,timesteps))
    dt = total_time/(timesteps)
    wdot_array = np.matrix(np.zeros((3,timesteps)))
    w_array = np.matrix(np.zeros((3,timesteps)))
    euler_array = np.matrix(np.zeros((3,timesteps)))
    orbital_quaternion = 1 #TODO: find what to do with this



    #configuration
    #initial euler angles
    euler_array[:,0]    = np.matrix([ [0], [0], [0] ])  #assumed initially in line with o frame (LVLH)
    #inital angular velocity
    w_array[:,0]        = np.matrix([ [0.01], [0.1], [0.01] ]) 

    #satellite physical characteristics
    Jmat       = np.matrix([     #kgm^2
        [ 1440, 0   , 0   ],
        [ 0   , 2500, 0   ],
        [ 0   , 0   , 3845]
    ])
    Jmat_inv    = np.linalg.inv(Jmat)
    wheel_h     = np.matrix([
        [0],
        [1000],
        [0]
    ]) #kgm^2/s


    #numerical solution
    for i in range(timesteps-1):
        #calculate n+1 wdot
        wdot_array[:,i+1]   = Get_Wdot(euler_array[:,i], w_array[:,i], Jmat, Jmat_inv, wheel_h, w_oi, u, R)
        #calculate n+1 w
        w_array[:,i+1]      = w_array[:,i] + np.multiply( dt, wdot_array[:,i+1])
        #update euler angles 
        euler_array[:,i+1]  = euler_array[:,i] + np.multiply( dt, Get_Euler_Rates(euler_array[:,i], w_array[:,i+1]) ).transpose()
        #repeat from here

    array   = np.concatenate(( euler_array, w_array, wdot_array ))
    df      = pd.DataFrame(array)
    df.to_excel(excel_writer = "/home/frosty/Projects/test.xlsx")

    plt.plot( t_vec, euler_array[0,:], 'r', t_vec, euler_array[1,:], 'g', t_vec, euler_array[2,:], 'b')
    plt.show()

    
