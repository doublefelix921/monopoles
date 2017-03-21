import math
import numpy as np

#X-field parameters
B_X=4.6e-10
Theta0_X=49.0
rc_X=4.8
r_X=2.9

def mag(V):
    return math.sqrt(V[0]**2+V[1]**2+V[2]**2)

pos=np.array([-8.5,8.661457517e-12,-9.72233787e-11])
r=math.sqrt(pos[0]**2+pos[1]**2)

def R_p(R,Z):
    if abs(Z) >= math.tan(Theta0_X)*(R-rc_X):
        return R*rc_X/(rc_X+abs(Z)/math.tan(Theta0_X))
    else:
        return R-abs(Z)/math.tan(Theta0_X)
        
r_p=R_p(r,pos[2])

bhat_X=(pos[2]**2+(r-r_p)**2)**(-1/2)*np.array([pos[0]*(r-r_p)/r, pos[1]*(r-r_p)/r, pos[2]])
print("z",pos[2])
print("r",r)
print("r_p",r_p)
print("z^2",pos[2]**2)
print("(r-r_p)**2",(r-r_p)**2)
print("sqrt(z^2+(r-r_p)^2",math.sqrt(pos[2]**2+(r-r_p)**2))
print(1/(pos[2]**2+(r-r_p)**2)**(-1/2))
print(mag(np.array([pos[0]*(r-r_p)/r, pos[1]*(r-r_p)/r, pos[2]])))
print(bhat_X, mag(bhat_X))