#This file will contain general functions about mathmatics of rotation
from matplotlib.pyplot import thetagrids
import numpy as np


def Angularvelocity2Q_dot(angular_velocity, quaternion):
    """
    I think this should be passive... what exactly does this mean?
    passive would make more sense... you want to describe the motion of the frame for attitude dynamics here.
    """
    q_dot = np.matmul( Q2PHI(quaternion), angular_velocity) 
    q_dot = np.multiply( 0.5, q_dot )
    return np.array([ q_dot[0,0], q_dot[0,1], q_dot[0,2], q_dot[0,3] ])

def Quat2Euler(quaternion):  #transform quaternions to 3-2-1 euler angles TODO support other euler angle rotations
    DCM = Quat2DCM(quaternion)
    phi = np.arctan2( DCM[1,2], DCM[2,2] )
    theta = np.sin( DCM[0,2] )
    psi = np.arctan2( DCM[0,1], DCM[0,0] )
    return phi, theta, psi

    
def Quat2DCM(quaternion):
    DCM = np.matmul( np.transpose(Q2PHI(quaternion)), Q2PSI(quaternion) )
    return DCM


def Q_Tran_Vec(quaternion, vector, active=False):
    DCM = np.matmul( np.transpose(Q2PHI(quaternion)), Q2PSI(quaternion) )

    if active==True:
        v_prime = np.matmul( np.transpose(DCM), vector ) 
    elif active==False:
        v_prime = np.matmul( DCM, vector ) #I THINK TRANSPOSE IS ACTIVE.. NO TRANSPOSE DCM IS PASSIVE

    return np.array([ v_prime[0,0], v_prime[0,1], v_prime[0,2] ])


def Q2PHI(quaternion):
    q = quaternion
    PHI = np.matrix([
        [  q[3], -q[2],  q[1]  ],
        [  q[2],  q[3], -q[0]  ],
        [ -q[1],  q[0],  q[3]  ],
        [ -q[0], -q[1], -q[2]  ],
    ])
    return PHI


def Q2PSI(quaternion):
    q = quaternion
    PSI = np.matrix([
        [  q[3],  q[2], -q[1]  ],
        [ -q[2],  q[3],  q[0]  ],
        [  q[1], -q[0],  q[3]  ],
        [ -q[0], -q[1], -q[2]  ],
    ])
    return PSI


if __name__=="__main__":
    
    #show identity quaternion and what it does to a vector
    q = np.array([ 0, 0 , 0, 1.0 ])  #identity quatenion i.e. no rotation of vector; equates to 1 as scalar part
    vector = np.array([ 1, 1, 1 ])
    v_prime = Q_Tran_Vec( q, vector )  #default to passive transformation
    print( v_prime )

    #show identity quaternion and what it does to a vector
    a = 1/np.sqrt(2)
    q = np.array([ 0, 0, a, a ])  #identity quatenion i.e. no rotation of vector; equates to 1 as scalar part
    vector = np.array([ 1, 1, 1 ])
    v_prime = Q_Tran_Vec( q, vector )  #default is passive transformation
    print( v_prime )

    #example actual quaternion rotation 
    q = np.array([ 0.043, -0.328, 0.764, 0.554 ])
    vector = np.array([ -0.945, 0.277, 0.172 ]) 
    v_prime = Q_Tran_Vec( q, vector, True)  #3RD True arg makes the transformation active.
    print( v_prime ) #answer should be [ 0.067, -0.915, -0.397]

    #testing attitude rate function
    q = np.array([ 0, 0, a, a])  #initial rotation of 90 deg about z axis
    w = np.array([ 1, 1, 1 ])  #angular velocity about current body axis... test value for now
    q_dot = Angularvelocity2Q_dot( w, q)
    print( q_dot )

