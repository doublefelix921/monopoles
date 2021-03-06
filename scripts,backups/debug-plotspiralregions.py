from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import math

fig = plt.figure()
ax = fig.gca(projection='3d')

step=1
minx=-20
maxx=20
miny=-20
maxy=20

x = np.arange(minx, maxx, step)
y = np.arange(miny, maxy, step)
x, y = np.meshgrid(x, y)



#=======CONSTANT PARAMETERS=======
c=9.715611713408621e-12 #light speed in kpc/s
exitstatus=True         #did it leave the 20kpc sphere by the end of the program?
nameofscript="final"    #name of the output file

#=======DEFAULT STARTING CONDITIONS=======
q_b=1.06655979e-28                                                   #the magnetic charge of the particle in Ampere*kpc by dirac quantization cond.
m=(10**13)*(1.78266184*10**-27)                                      #mass (kg, enter in GeV). An estimate from Wick (2002).
pos_0=np.array([-8.5, 0.0, 0.0])                                     #starting position (kpc), position of earth is x=-8.33937979kpc
#vel_0=np.array([0.0, 0.0, -300])*(3.24077929e-20)                   #starting velocity (kpc/s), but enter the numbers in m, conversion is there.
dt=1.0e11                                                            #the timestep (smaller = more accuracy, more computing time) in seconds
distance_to_halt= 50.0                                               #when to stop the program if it doesn't leave r>20kpc
V=300*(3.24077929e-20)
posx=pos_0[0]
posy=pos_0[1]
posz=pos_0[2]
THETA=0.0
PHI=0.0

#=======DISK PARAMETERS, WITH B IN TESLA=======
b1=0.1e-10 
b2=3.0e-10
b3=-0.9e-10
b4=-0.8e-10
b5=-2.0e-10
b6=-4.2e-10
b7=0.0e-10
b8=2.7e-10
b_ring=0.1e-10
h_disk=0.4
w_disk=0.27

#=======HALO PARAMETERS & FUNCTIONS=======
B_n=1.4e-10
B_s=-1.1e-10
r_n=9.22
r_s=16.7 #NOTE: parameter has a large error.
w_h=0.2
z_0=5.3
i=math.radians(11.5)                                          #this is the opening "pitch" angle of the logarithmic spiral boundaries.
r_negx=np.array([5.1, 6.3, 7.1, 8.3, 9.8, 11.4, 12.7, 15.5])  #these are the radii at which the spirals cross the x-axis

def r_i(T,I):                                                 #halo spiral boundary function
    return r_negx[I]*math.e**((T-math.pi)*math.tan(i))
    
def L(Z,H,W):                                                 #disk-halo transition function
    return (1.0+math.e**(-2.0*(abs(Z)-H)/W))**-1.0

#=======X-FIELD PARAMATERS & FUNCTIONS=======
B_X=4.6e-10
Theta0_X=math.radians(49.0)
rc_X=4.8
r_X=2.9

def R_p(R,Z):
    if abs(Z) >= math.tan(Theta0_X)*(R-rc_X):
        return R*rc_X/(rc_X+abs(Z)/math.tan(Theta0_X))

    else:
        return R-abs(Z)/math.tan(Theta0_X)

def b_Xf(R_P):
    return B_X*math.e**(-R_P/r_X)

#=======OTHER PRELIMINARY DEFINITIONS=======

def mag(V):
    return math.sqrt(V[0]**2+V[1]**2+V[2]**2)
    
def Theta(X,Y):
    if X == 0.0:
        if Y > 0.0:
            return math.pi/2
        if Y < 0.0:
            return 3.0*math.pi/2
    else:
        if Y < 0.0:
            return math.atan2(Y,X)+2*math.pi
        else:
            return math.atan2(Y,X)
                      
def region(pos):
    r=math.sqrt(pos[0]**2.0+pos[1]**2.0)
    theta=Theta(pos[0],pos[1])
    r_p=R_p(r,pos[2])
    b_X=b_Xf(r_p)
    
    
    bfield=np.array([0.0, 0.0, 0.0]) #0 for r>1, <20
    if mag(pos) > 1.0:
        #halo:
#        if pos[2] >= 0.0:
#            bfield = math.e**(-abs(pos[2])/z_0)*L(pos[2],h_disk,w_disk)* B_n*(1-L(r,r_n,w_h)) * np.array([-1*pos[1]/r, pos[0]/r, 0.0])
#        else:
#            bfield = math.e**(-abs(pos[2])/z_0)*L(pos[2],h_disk,w_disk)* B_s*(1-L(r,r_s,w_h)) * np.array([-1*pos[1]/r, pos[0]/r, 0.0])
            
#        #X-field: 
#        if pos[2] != 0: 
#            bhat_X=(pos[2]**2+(r-r_p)**2)**(-0.5)*np.array([pos[0]*(r-r_p)/r, pos[1]*(r-r_p)/r, pos[2]])
#            #if pos[2] < 0.0:
#            #    bhat_X[2]=-1.0*bhat_X[2]
#            if abs(r_p) <  rc_X:
#                bfield += b_X*(r_p/r)**2*bhat_X
#            if abs(r_p) >= rc_X:
#                bfield += b_X*(r_p/r)*bhat_X
#                
        #disk:
        if r>= 3.0 and r < 5.0 :
            bfield+=b_ring*(1-L(pos[2],h_disk,w_disk))
            return -1.0
            
        elif r>=5.0 and r<=20.0:
            if   r >= r_i(theta,7) and r < r_i(theta,0): #region 1
                bfield+=(b1/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
                print("r1")
                return 1.0
            elif r >= r_i(theta,0) and r < r_i(theta,1): #region 2
                bfield+=(b2/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
                print("r2")
                return 2.0
            elif r >= r_i(theta,1) and r < r_i(theta,2): #region 3
                bfield+=(b3/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
                print("r3")
                return 3.0
            elif r >= r_i(theta,2) and r < r_i(theta,3): #region 4
                bfield+=(b4/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
                print("r4")
                return 4.0
            elif r >= r_i(theta,3) and r < r_i(theta,4): #region 5
                bfield+=(b5/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
                print("r5")
                return 5.0
            elif r >= r_i(theta,4) and r < r_i(theta,5): #region 6
                bfield+=(b6/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
                print("r6")
                return 6.0
            elif r >= r_i(theta,5) and r < r_i(theta,6): #region 7
                bfield+=(b7/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
                print((b7/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0]))
                print "b7",b7,"\nz",pos[2],"\nh_disk",h_disk,"\nw_disk",w_disk,"\nr",r,"\ntheta",theta,"\ni",i,"\nL",L(pos[2],h_disk,w_disk)
                print
                print "direction",np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
                print("region r7")
                print"bfield",bfield
                return 7.0
            elif r >= r_i(theta,6) and r < r_i(theta,7): #region 8
                bfield+=(b8/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
                print("r7")
                return 8.0



#print(calculate_bfield([-5,-5,0]))
#R = np.sqrt(X**2 + Y**2)
#Z = np.sin(R)
Rvalues=np.zeros(shape=(len(x),len(y)))
z=0.0

for i in range(len(x)):
    for j in range(len(y)):
        POSITION=np.array([i*step+minx, j*step+miny, z])
        Rvalues[i][j]=region(POSITION)
        #print(POSITION)
        #print(Bvalues[i][j])
        #print(mag(calculate_bfield(pos)))

surf = ax.plot_surface(x, y, Rvalues, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
#ax.set_zlim(-1.01, 1.01)

#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#for i in range(len(Bvalues)):
#    print(Bvalues[i])

fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()

