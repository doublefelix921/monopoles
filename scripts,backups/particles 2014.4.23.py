#Tracking Magnetic Monopoles Through the Galaxy:
#input starting conditions into the program and it will track the magnetic monopole's path through the bfield.

#Bfield given by "A New Model of the Galatic Magnetic Field" Jansson, Farrar (2012)
#Mass bounds from "Signatures for a Cosmic Flux of Magnetic Monopoles" Wick (2002) are 40TeV <~ M <~ 10^8TeV
#Wick(2002) has more information on likely mass values

#units used are all SI units except: distance is in kpc, and angles are in degrees, magnetic field strength is in microgauss.

#assumptions made: 
#1.the magnetic monopoles have no electric charge.
#2.there are no collisions - space is empty.
#3.the electric and gravitational fields everywhere are negligible. No forces are involved except the magnetic force.

import math
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys
import random

#starting conditions:
q_b=5.19230065e-39                                                      #the magnetic charge of the particle in Ampere*kpc by dirac quantization cond.
m=1.78266184e-10                                                        #mass (kg). An estimate from Wick (2002).
pos_0=np.array([-8.33937979, 0.00, 0.00])                               #starting position (kpc)
vel_0=np.array([0.0, 0.0, 0.0])*(3.24077929e-20)                       #starting velocity (kpc/s), but enter the numbers in m, conversion is there.
dt=1e8                                                                 #the timestep (smaller = more accuracy, more computing time) in seconds
distance_to_track= .00001                                                    #how far to track the particle (kpc)

#disk parameters. B is in tesla.
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

#halo parameters
B_n=1.4e-10
B_s=-1.1e-10
r_n=9.22
r_s=16.7 #NOTE: parameter has a large error.
w_h=0.2
z_0=5.3

#X-field parameters
B_X=4.6e-10
Theta0_X=49.0
rc_X=4.8
r_X=2.9

#striation parameters (the striated field is currently unused)
#gamma= 2.92
#alpha= 2.65 #taken from page 9, top right paragraph
#beta= gamma/alpha-1

#other preliminary measures:
c=9.715611713408621e-12    #speed of light in kpc/s

def mag(V):
    return math.sqrt(V[0]**2+V[1]**2+V[2]**2)
    
def Theta(X,Y):
    if X=0:
        if Y>0:
            return math.pi/2
        if Y<0:
            return 3*math.pi/2
    else:
        return math.atan(Y/X)
        
pos=np.array([0.0, 0.0, 0.0])
for N in range(2):
    pos[N]=pos_0[N]
vel=vel_0

trailx=(pos[0],)                    #trailx,y,z is used to save the coordinates of the particle at each step, to plot the path afterwards
traily=(pos[1],)
trailz=(pos[2],)

gam=1/math.sqrt(1-(mag(vel)/c)**2)

if gam != 1:
    KE=m*c**2*(gam-1)*5.942795e48       #KE, converted to GeV
else:
    KE=m*mag(vel)**2/2
KEhistory=(KE,)

distance_tracked=0.0                  #set the distance travelled so far to 0

time=0.0

#boundary function (between disk and halo fields)
def L(Z,H,W):
    return (1+math.e**(-2*(abs(Z)-H)/W))**-1

#halo boundary spirals:
i=11.5                                                        #this is the opening "pitch" angle of the logarithmic spiral boundaries.
r_negx=np.array([5.1, 6.3, 7.1, 8.3, 9.8, 11.4, 12.7, 15.5])  #these are the radii at which the spirals cross the x-axis

def r_i(T,I):
    return r_negx[I]*math.e**((T-math.pi)*math.tan(math.radians(i)))
    
#X-field definitions:

def r_p(R,Z):
    if abs(Z) >= math.tan(Theta0_X)*(R-rc_X):
        return R*rc_X/(rc_X+abs(Z)/math.tan(Theta0_X))
    else:
        return R-abs(Z)/math.tan(Theta0_X)

def b_X(R_P):
    return B_X*math.e**(-R_P/r_X)

#def Theta_X(R,Z):
#    if abs(Z)>=math.tan(Theta0_X)*math.sqrt(R)-math.tan(Theta0_X)*rc_X:
#        return math.atan(abs(Z)/(R-r_p(R,Z)))
#    else:
#        return Theta0_X
    
#preliminary check:
if mag(vel) >= c:
    print("Error: Initial velocity cannot exceed the speed of light. Currently, it is",mag(vel)/c,"times c.")
    sys.exit("Stopped program.")
    

#print initial conditions:
print()
print()
print("=========================PARTICLE TRAIL INFO=========================")
print("Your Initial Parameters: \nq_b =",q_b,"A*kpc     m =",m,"kg     dt =",dt,"s     distance to track =",distance_to_track,"kpc","KE =",KE,"GeV")
print("initial position (kpc)   =",pos_0,"\ninitial velocity (kpc/s) =",vel_0)
print()



#ok, let's start tracking the monopole. Most of the calculations in the loop are to find the bfield.

while distance_tracked < distance_to_track:
    #definitions:
    r=math.sqrt(pos[0]**2+pos[1]**2)
    theta=Theta(pos[0],pos[1])
    gam=1/math.sqrt(1-(mag(vel)/c)**2)
    
    #now for bfield calculation, component by component: halo, then X-field, then disk (striated not currently used)
        
    #halo component:
    if pos[2]>=0:
        bfield = math.e**(-abs(pos[2])/z_0)*L(pos[2],h_disk,w_disk)* B_n*(1-L(r,r_n,w_h)) * np.array([-1*pos[1]/r, pos[0]/r,0])
    else:
        bfield = math.e**(-abs(pos[2])/z_0)*L(pos[2],h_disk,w_disk)* B_s*(1-L(r,r_s,w_h)) * np.array([-1*pos[1]/r, pos[0]/r,0])
        
    #X-field component:
    
    if r_p(r,pos[2]) < rc_X:
        bfield += b_X(r_p(r,pos[2]))*(r_p(r,pos[2])/r)**2*((r-r_p(r,pos[2]))**2+pos[2]**2)**(-1/2)*np.array([(1-r_p(r,pos[2])/r)*pos[0],(1-r_p(r,pos[2])/r)*pos[1],pos[2]])
    else:
        bfield += b_X(r_p(r,pos[2]))*(r_p(r,pos[2])/r)*((r-r_p(r,pos[2]))**2+pos[2]**2)**(-1/2)*np.array([(1-r_p(r,pos[2])/r)*pos[0],(1-r_p(r,pos[2])/r)*pos[1],pos[2]])

    #disk component:
    if r>= 3.0 and r < 5.0 :
        bfield+=b_ring*(1-L(pos[2],h_disk,w_disk))
        
    elif r>=5.0 and r<=20.0:
        if   r >= r_i(theta,7) and r < r_i(theta,0): #region 1
            bfield+=(b1/r)*(1-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+math.radians(i)),-math.cos(theta+math.radians(i)),0])
            
        elif r >= r_i(theta,0) and r < r_i(theta,1): #region 2
            bfield+=(b2/r)*(1-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+math.radians(i)),-math.cos(theta+math.radians(i)),0])
            
        elif r >= r_i(theta,1) and r < r_i(theta,2): #region 3
            bfield+=(b3/r)*(1-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+math.radians(i)),-math.cos(theta+math.radians(i)),0])
            
        elif r >= r_i(theta,2) and r < r_i(theta,3): #region 4
            bfield+=(b4/r)*(1-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+math.radians(i)),-math.cos(theta+math.radians(i)),0])
            
        elif r >= r_i(theta,3) and r < r_i(theta,4): #region 5
            bfield+=(b5/r)*(1-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+math.radians(i)),-math.cos(theta+math.radians(i)),0])
            
        elif r >= r_i(theta,4) and r < r_i(theta,5): #region 6
            bfield+=(b6/r)*(1-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+math.radians(i)),-math.cos(theta+math.radians(i)),0])
            
        elif r >= r_i(theta,5) and r < r_i(theta,6): #region 7
            bfield+=(b7/r)*(1-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+math.radians(i)),-math.cos(theta+math.radians(i)),0])
            
        elif r >= r_i(theta,6) and r < r_i(theta,7): #region 8
            bfield+=(b8/r)*(1-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+math.radians(i)),-math.cos(theta+math.radians(i)),0])
                                                                    
    if r<1 or r>20:
        bfield=0*bfield #not enough data to extrapolate field, it was set to 0 near galactic center or past r=20kpc.
        
    #striated fields (unfinished, unused):
    #zeroorone=randrange(2)
    #if zeroorone==0:
    #    bfield-= math.sqrt(beta)*bfield
    #if zeroorone==1:
    #    bfield+= math.sqrt(beta)*bfield
    
    #CALCULATION OF THE CHANGE IN POSITION:
    #nonrelativistic:
    #acc=bfield*(q_b/m)
    
    #pos=np.array([pos[0]-vel[0]*dt-0.5*acc[0]*(dt**2),
    #              pos[1]-vel[1]*dt-0.5*acc[1]*(dt**2),
    #              pos[2]-vel[2]*dt-0.5*acc[2]*(dt**2)])
    #distance_tracked+=math.sqrt((vel[0]*dt+0.5*acc[0]*(dt**2))**2+(vel[1]*dt+0.5*acc[1]*(dt**2))**2+(vel[2]*dt+0.5*acc[2]*(dt**2))**2)
    
    #vel=np.array([vel[0]-acc[0]*dt,
    #              vel[1]-acc[1]*dt,
    #              vel[2]-acc[2]*dt])
    
    #trailx=trailx+(pos[0],)
    #traily=traily+(pos[1],)
    #trailz=trailz+(pos[2],)
    
    
    
    #KE=9.521406e38*6.24150934e15*gam*m*c**2 #calculate KE, and convert from kg*kpc^2/s^2 to J to MeV
    #KEhistory=KEhistory+(KE,)
    
    #RELATIVISTIC:
        
    force=q_b*bfield #In kg*kpc/s^2
    
    acc_prefactor=(gam*m*(c**2+gam**2*mag(vel)**2))**-1
    acc=np.array([acc_prefactor*(force[0]*(c**2+gam**2*vel[1]**2+gam**2*vel[2]**2) - gam**2*vel[0]*(force[1]*vel[1]+force[2]*vel[2])),
                  acc_prefactor*(force[1]*(c**2+gam**2*vel[0]**2+gam**2*vel[2]**2) - gam**2*vel[1]*(force[0]*vel[0]+force[2]*vel[2])),
                  acc_prefactor*(force[2]*(c**2+gam**2*vel[0]**2+gam**2*vel[1]**2) - gam**2*vel[2]*(force[0]*vel[0]+force[1]*vel[1]))])
    
    vel_i=vel
    
    vel+= -acc*dt
    
    pos+= -vel_i*dt-0.5*acc*dt**2
    
    if gam!=1:
        KE=m*c**2*(gam-1)*5.942795e48       #KE, converted to GeV
    else:
        KE=m*mag(vel)**2/2
    KEhistory+=(KE,)
    
    time+=dt
    
    if random.randint(1,10000)==1:
        print("distance traveled:",distance_tracked)
        print("pos",pos)
        print("vel",vel)
        print("KE",KE)
        print("time",time)
    
    distance_tracked+=math.sqrt((vel_i[0]*dt+0.5*acc[0]*dt**2)**2+(vel_i[1]*dt+0.5*acc[1]*dt**2)**2+(vel_i[2]*dt+0.5*acc[2]*dt**2)**2)
    
    trailx=trailx+(pos[0],)
    traily=traily+(pos[1],)
    trailz=trailz+(pos[2],)
    

distance_from_start=math.sqrt( (pos[0]-pos_0[0])**2 +(pos[1]-pos_0[1])**2 +(pos[2]-pos_0[2])**2)

print("The final position (kpc) is ( " + str(pos[0]) + ", " + str(pos[1]) + ", " + str(pos[2]) + ")." )
print("The final velocity (kpc/s) is ( " + str(vel[0]) + ", " + str(vel[1]) + ", " + str(vel[2]) + ")." )
print("The final Kinetic Energy (GeV) is",KE)
print()
print("Distance from initial position is", distance_from_start,"kpc")
print("The journey took",time,"seconds.")
print()
print("The galactic center is plotted as a blue dot, and the sun is plotted as a yellow dot.")
print()
print()

fig = plt.figure(figsize=(5,8))
fig.canvas.set_window_title('Monopole Trail')

ax = fig.add_subplot(211, projection='3d')
ax.plot(trailx,traily,trailz)
ax.plot([-8.33937979],[0],[0],'yo')
#ax.plot([0],[0],[0],'bo')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.suptitle("Particle Trail (click and drag to change view)",fontsize=12,position=(0.5,0.93),weight='bold')
plt.title('Kinetic Energy (GeV)', position=(0.5,-0.2),fontsize=12,weight='bold')

t_array=np.arange(0,dt*len(KEhistory),dt)
ax2=fig.add_subplot(212)
ax2.plot(t_array,KEhistory,)
ax2.set_xlabel("Time (s)")
ax2.set_ylabel("Particle Kinetic Energy (GeV)")
plt.grid(True)

plt.show()
