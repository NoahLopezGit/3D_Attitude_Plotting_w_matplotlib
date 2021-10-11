#This file will contain the code used in HW 2
import numpy as np
import matplotlib.pyplot as plt
import attitudematrices as att
import THREED_Att_Animation 


#Problem 1 comparing analytical and numerical solutions of EOMS

#calculate one update step using eulers method
def Euler_Step(Yn, Ynprime, h): #input  Yn state and Yn' value, step size
    Ynp1 = Yn + h*Ynprime
    #return updated step
    return Ynp1 

def Update_Vec(vec, quat):
    v_tmp = att.Q_Tran_Vec( quat, vec )
    v_tmp = np.multiply( 1/np.linalg.norm(v_tmp), v_tmp) #normalizing vector bc this blows up TODO: find why this is unstable
    return v_tmp


if __name__=="__main__":
    
    #step size and timesteps
    h = 0.01
    timesteps = 1000
    t = np.linspace(0, h*(timesteps), timesteps+1)


    #defining derivatives as functions 
    w1_prime = lambda w2, w3 :  3/10*w2*w3
    w2_prime = lambda w1, w3 : -3/10*w1*w3
    w3_prime = 0

    #analytical solution
    w3_anal = 1
    w1_anal = lambda t : 0.1*np.cos(3/10*t) + 0.1*np.sin(3/10*w3_anal*t)
    w2_anal = lambda t : 0.1*np.cos(3/10*t) - 0.1*np.sin(3/10*w3_anal*t)
    
    #create list for angular velocties
    w1_est = [0.1] #first index contains initial state
    w2_est = [0.1]
    w3_est = [1]
    
    w1_exact = [0.1]
    w2_exact = [0.1]
    w3_exact = [1]

    for i in range(timesteps):
        w1_tmp = Euler_Step( w1_est[-1], w1_prime( w2_est[-1], w3_est[-1]), h )
        w2_tmp = Euler_Step( w2_est[-1], w2_prime( w1_est[-1], w3_est[-1]), h )
        w3_tmp = Euler_Step( w3_est[-1], w3_prime , h )

        w1_est.append( w1_tmp )
        w2_est.append( w2_tmp )
        w3_est.append( w3_tmp )

        w1_exact.append( w1_anal((i+1)*h) )
        w2_exact.append( w2_anal((i+1)*h) )
        w3_exact.append( w3_anal )  

    
    #solve for angular state as a function of time???

    quaternion = [np.array([ 0, 0, 0, 1 ])]  #initial identity quaternion
    vector_z = [np.array([ 0, 0, 1])]  #initial vector list
    vector_y = [np.array([ 0, 1, 0])] 
    vector_x = [np.array([ 1, 0, 0])]
    
    #propogation of attitude state
    for i in range(timesteps):

        #calculate new quaternion
        angular_velocity = np.array([ w1_est[i], w2_est[i], w3_est[i] ])
        q_dot = att.Angularvelocity2Q_dot( angular_velocity,  quaternion[-1])
        q_update = np.array([   Euler_Step( quaternion[-1][0], q_dot[0], h),  #TODO: I beleive this is causing the instability...find out why
                                Euler_Step( quaternion[-1][1], q_dot[1], h),
                                Euler_Step( quaternion[-1][2], q_dot[2], h),
                                Euler_Step( quaternion[-1][3], q_dot[3], h) ])
        q_update = np.multiply( 1/np.linalg.norm(q_update), q_update)  #normalizing quaternion
        quaternion.append( q_update )

        #calculate new body vectors with new quaternion
        vector_z.append(Update_Vec( vector_z[0], quaternion[-1]))  
        vector_y.append(Update_Vec( vector_y[0], quaternion[-1]))
        vector_x.append(Update_Vec( vector_x[0], quaternion[-1]))  

    
    #adjusting format of quaternion and vector lists
    quaternion_plot     = np.zeros( (4, timesteps+1) )
    phi                 = np.zeros( timesteps+1 ) 
    theta               = np.zeros( timesteps+1 ) 
    psi                 = np.zeros( timesteps+1 )
    
    data_z = np.zeros( (3, timesteps+1) )
    data_y = np.zeros( (3, timesteps+1) )
    data_x = np.zeros( (3, timesteps+1) )

    for i in range(timesteps+1):
        data_z[:,i]     = vector_z[i][:]
        data_y[:,i]     = vector_y[i][:]
        data_x[:,i]     = vector_x[i][:]
        quaternion_plot[:, i] = quaternion[i][:]
        #euler angles
        phi[i], theta[i], psi[i] = att.Quat2Euler(quaternion[i][:])


    fig, axs = plt.subplots(2, 2)
    #show angular velocity, quaternion, euler angles
    axs[0,0].plot(  t, w1_est, 'r', t, w2_est, 'g', t, w3_est, 'b',  
                    t, w1_exact, 'r--', t, w2_exact, 'g--', t, w3_exact, 'b--' )
    axs[0,1].plot(  t, quaternion_plot[0,:], 'r', t, quaternion_plot[1,:], 'g', 
                    t, quaternion_plot[2,:], 'b', t, quaternion_plot[3,:], 'y' )
    axs[1,1].plot( t, phi, 'r', t, theta, 'g', t, psi, 'b')

    plt.show()
    
    #sending data to 3d animator; this will play in 0.1x real speed 
    THREED_Att_Animation.Animate_Attitude_Set(data_y, data_x, data_z, h)
