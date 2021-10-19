#Example for attitude kinematics module
import Attitude_Kinematics as att
import numpy as np
import matplotlib.pyplot as plt



#step size and timestepsS
h = 0.25
timesteps = 40
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
    w1_tmp = att.Euler_Step( w1_est[-1], w1_prime( w2_est[-1], w3_est[-1]), h )
    w2_tmp = att.Euler_Step( w2_est[-1], w2_prime( w1_est[-1], w3_est[-1]), h )
    w3_tmp = att.Euler_Step( w3_est[-1], w3_prime , h )

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
    q_update = att.Update_Q( quaternion[i], angular_velocity, h )
    q_update = np.multiply( 1/np.linalg.norm(q_update), q_update)  #normalizing quaternion
    quaternion.append( q_update )

    #calculate new body vectors with new quaternion
    vector_z.append(att.Update_Vec( vector_z[0], quaternion[-1]))  
    vector_y.append(att.Update_Vec( vector_y[0], quaternion[-1]))
    vector_x.append(att.Update_Vec( vector_x[0], quaternion[-1]))  


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
att.Animate_Attitude_Set(data_y, data_x, data_z, h)


plt.plot(   t, w1_est, 'r', t, w2_est, 'g', t, w3_est, 'b',  
            t, w1_exact, 'r--', t, w2_exact, 'g--', t, w3_exact, 'b--' )
plt.title( "Exact Angular Rates (dashed) vs Estimated (solid)" )
plt.xlabel("Time (s)")
plt.show()