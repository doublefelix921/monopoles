from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import math
import datetime

fig = plt.figure()
ax = fig.gca(projection='3d')

x, y, z = np.meshgrid(np.arange(-20, 20, 1),
                      np.arange(-20, 20, 1),
                      np.arange(-20, 20, 1))


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
posx=x
posy=y
posz=z
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

#def Theta_X(R,Z):
#    if abs(Z)>=math.tan(Theta0_X)*math.sqrt(R)-math.tan(Theta0_X)*rc_X:
#        return math.atan(abs(Z)/(R-r_p(R,Z)))  #think about atan and maybe problems with quadrants
#    else:
#        return Theta0_X

#striation parameters (the striated field is currently unused)
#gamma= 2.92
#alpha= 2.65 #taken from page 9, top right paragraph
#beta= gamma/alpha-1

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
    
def Gam(V):
    return (1.0-(mag(V)/c)**2)**(-0.5)
    
    #check the KE
def KE(VEL):
    if Gam(VEL) != 1:
        KE=m*c**2*(Gam(VEL)-1)*5.942795e48     #KE, converted to GeV
    else:
        KE=m*mag(VEL)**2*5.942795e48/2
        

def gettimestd():
    output=""
    W=str(datetime.datetime.now())
    for n in [2,3,5,6,8,9]:
        output+=W[n]
    output+="_"
    for n in range(11,13):
        output+=W[n]
    output+="."
    for n in range(14,16):
        output+=W[n]
    output+="."
    for n in range(17,23):
        output+=W[n]
    return output
    
def gettimesecs():
    W=str(datetime.datetime.now())
    tsecs=0.0
    tsecs+=float(W[0]+W[1]+W[2]+W[3])*31556926
    tsecs+=float(W[5]+W[6])*2629743.83
    tsecs+=float(W[8]+W[9])*86400
    tsecs+=float(W[11]+W[12])*3600
    tsecs+=float(W[14]+W[15])*60
    tsecs+=float(W[17]+W[18])*1
    tsecs+=float("0."+W[20]+W[21]+W[22]+W[23]+W[24])
    return tsecs
    
tzero=gettimesecs()

g=open(nameofscript+"_full_output"+str(gettimestd())+".txt","a")

bfield=np.array([0.0, 0.0, 0.0]) #0 for r>1, <20
if math.sqrt(x**2+y**2+z**2) > 1.0:
    #halo:
    if z >= 0.0:
        bfield = math.e**(-abs(z)/z_0)*L(z,h_disk,w_disk)* B_n*(1-L(r,r_n,w_h)) * np.array([-1*y/r, x/r, 0.0])
    else:
        bfield = math.e**(-abs(z)/z_0)*L(z,h_disk,w_disk)* B_s*(1-L(r,r_s,w_h)) * np.array([-1*y/r, x/r, 0.0])
        
    #X-field: 
    if z != 0: 
        bhat_X=(z**2+(r-r_p)**2)**(-0.5)*np.array([x*(r-r_p)/r, y*(r-r_p)/r, z])
        if z < 0.0:
            bhat_X[2]=-1.0*bhat_X[2]
        if abs(r_p) <  rc_X:
            bfield += b_X*(r_p/r)**2*bhat_X
        if abs(r_p) >= rc_X:
            bfield += b_X*(r_p/r)*bhat_X
            
    #disk:
    if r>= 3.0 and r < 5.0 :
        bfield+=b_ring*(1-L(z,h_disk,w_disk))
        
    elif r>=5.0 and r<=20.0:
        if   r >= r_i(theta,7) and r < r_i(theta,0): #region 1
            bfield+=(b1/r)*(1.0-L(z,h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
            
        elif r >= r_i(theta,0) and r < r_i(theta,1): #region 2
            bfield+=(b2/r)*(1.0-L(z,h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
            
        elif r >= r_i(theta,1) and r < r_i(theta,2): #region 3
            bfield+=(b3/r)*(1.0-L(z,h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
            
        elif r >= r_i(theta,2) and r < r_i(theta,3): #region 4
            bfield+=(b4/r)*(1.0-L(z,h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
            
        elif r >= r_i(theta,3) and r < r_i(theta,4): #region 5
            bfield+=(b5/r)*(1.0-L(z,h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
            
        elif r >= r_i(theta,4) and r < r_i(theta,5): #region 6
            bfield+=(b6/r)*(1.0-L(z,h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
            
        elif r >= r_i(theta,5) and r < r_i(theta,6): #region 7
            bfield+=(b7/r)*(1.0-L(z,h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
            
        elif r >= r_i(theta,6) and r < r_i(theta,7): #region 8
            bfield+=(b8/r)*(1.0-L(z,h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
            







u = np.sin(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)
v = -np.cos(np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)
w = (np.sqrt(2.0 / 3.0) * np.cos(np.pi * x) * np.cos(np.pi * y) *
     np.sin(np.pi * z))





ax.quiver(x, y, z, u, v, w, length=0.1)

plt.show()