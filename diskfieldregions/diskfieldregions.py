import matplotlib.pyplot as plt
import numpy as np
import math

theta=np.arange(-2*math.pi,4*math.pi,0.01)
five = np.ones_like(theta)*5
twenty = np.ones_like(theta)*20

ax1=plt.subplot(111,polar=True)
ax1.set_rlim((0,20))

i=math.radians(11.5)                                          #this is the opening "pitch" angle of the logarithmic spiral boundaries.
r_negx=np.array([5.1, 6.3, 7.1, 8.3, 9.8, 11.4, 12.7, 15.5])  #these are the radii at which the spirals cross the x-axis

def r_i(T,I):                                                 #halo spiral boundary function
    return r_negx[I]*math.e**((T-math.pi)*math.tan(i))
    
r1,r2,r3,r4,r5,r6,r7,r8=[[] for _ in range(8)]
for n in range(len(theta)):
    r1+=[r_i(theta[n],0)]
    r2+=[r_i(theta[n],1)]
    r3+=[r_i(theta[n],2)]
    r4+=[r_i(theta[n],3)]
    r5+=[r_i(theta[n],4)]
    r6+=[r_i(theta[n],5)]
    r7+=[r_i(theta[n],6)]
    r8+=[r_i(theta[n],7)]

ax1.plot(theta,r1,'b-')
ax1.plot(theta,r2,'g-')
ax1.plot(theta,r3,'r-')
ax1.plot(theta,r4,'k-')
ax1.plot(theta,r5,'b--')
ax1.plot(theta,r6,'g--')
ax1.plot(theta,r7,'r--')
ax1.plot(theta,r8,'k--')
ax1.plot(theta,five,linewidth=4)
ax1.plot(theta,twenty,linewidth=4)
#ax1.plot(t,math.cos(t),linewidth=4)

plt.show()
