import matplotlib.pyplot as plt
import numpy as np
import math

t=np.arange(-2*math.pi,8*math.pi,0.01)
five = np.ones_like(t)*5
twenty = np.ones_like(t)*20

ax1=plt.subplot(111,polar=True)
ax1.set_rlim((0,20))

arg=(t-math.pi)*math.tan(math.radians(11.5))

ax1.plot(t,5.1*math.e**(arg),'b-')
ax1.plot(t,6.3*math.e**(arg),'g-')
ax1.plot(t,7.1*math.e**(arg),'r-')
ax1.plot(t,8.3*math.e**(arg),'k-')
ax1.plot(t,9.8*math.e**(arg),'b--')
ax1.plot(t,11.4*math.e**(arg),'g--')
ax1.plot(t,12.7*math.e**(arg),'r--')
ax1.plot(t,15.5*math.e**(arg),'k--')
ax1.plot(t,five,linewidth=4)
ax1.plot(t,twenty,linewidth=4)
#ax1.plot(t,math.cos(t),linewidth=4)

plt.show()
