#Tracking Magnetic Monopoles Through the Galaxy:
#input starting conditions into the program and it will track the magnetic monopole's path through the bfield.

#Bfield given by "A New Model of the Galatic Magnetic Field" Jansson, Farrar (2012).
#Mass bounds from "Signatures for a Cosmic Flux of Magnetic Monopoles" Wick (2002) are 40TeV <~ M <~ 10^8TeV.
#Wick(2002) has more information on likely mass values.

#units used are all SI units except: distance is in kpc, and angles are in degrees, magnetic field strength is in microgauss.

#variable names generally match those in the Bfield paper.

from __future__ import print_function, division
import math
import numpy as np
import sys
import random
import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#=======CONSTANT PARAMETERS=======
c=9.715611713408621e-12 #light speed in kpc/s
exitstatus=True         #did it leave the 20kpc sphere by the end of the program?
nameofscript="mpoles15_graphed"    #name of the output file

microGtoT=1e-10
GeVtokg=1.78266184e-27
mtokpc=3.24077929e-20

#=======DEFAULT STARTING CONDITIONS=======
q_b=1.06655979e-28                                                  #the magnetic charge of the particle in Ampere*kpc by dirac quantization cond.
m=1e13*GeVtokg                                                      #mass (kg, enter in GeV). An estimate from Wick (2002).
dt=1.0e19                                                           #the timestep (smaller = more accuracy, more computing time) in seconds
distance_to_halt= 40.0                                              #when to stop the program if it doesn't leave r>20kpc

pos_0=np.array([-8.5, 0.0, 0.0])                                    #starting position (kpc), position of earth is x=-8.33937979kpc
#vel_0=np.array([0.0, 0.0, -300])*mtokpc                            #starting velocity (kpc/s), but enter the numbers in m, conversion is there.
V=0.0*mtokpc
THETA=0.0    #theta is the polar angle
PHI=0.0      #phi is the azimuthal angle



#=======INPUTTED STARTING CONDITIONS=======
#if sys.argv[1]!="N":
#    q_b=float(sys.argv[1])
#if sys.argv[2]!="N":
#    m=float(sys.argv[2])
#if sys.argv[3]!="N":
#    posx=float(sys.argv[3])
#if sys.argv[4]!="N":
#    posy=float(sys.argv[4])
#if sys.argv[5]!="N":
#    posz=float(sys.argv[5])
#if sys.argv[6]!="N":
#    distance_to_halt=float(sys.argv[6])
#if sys.argv[7]!="N":
#    dt=float(sys.argv[7])
#if sys.argv[8]!="N":
#    V=float(sys.argv[8])
#if sys.argv[9]!="N":
#    THETA=float(sys.argv[9])
#if sys.argv[10]!="N":
#    PHI=float(sys.argv[10])

#=======DISK PARAMETERS, WITH B IN microGauss. Converted to Tesla when field is calculated=======
b1=0.1
b2=3.0
b3=-0.9
b4=-0.8
b5=-2.0
b6=-4.2
b7=0.0
b8=2.7
b_ring=0.1
h_disk=0.4
w_disk=0.27

#=======HALO PARAMETERS & FUNCTIONS=======
B_n=1.4
B_s=-1.1
r_n=9.22
r_s=16.7 #NOTE: parameter has a large error.
w_h=0.2
z_0=5.3
i=np.radians(11.5)                                          #this is the opening "pitch" angle of the logarithmic spiral boundaries.
r_negx=np.array([5.1, 6.3, 7.1, 8.3, 9.8, 11.4, 12.7, 15.5])  #these are the radii at which the spirals cross the x-axis

def r_i(T,I):                                                 #halo spiral boundary function
    return r_negx[I]*np.e**((T-np.pi)*np.tan(i))
    
def subregioncheck(r,theta,i1,i2):
    if r>r_i(theta,i1) and r<r_i(theta,i2):
        return True
    else:
        return False

#the following function answers: am I in between regions i1 and i2?
def regioncheck(r,theta,i1,i2): #checks if the particle is in this region.
    if i1 != 7:
        if r>r_i(theta-2*np.pi,i1) and r<r_i(theta-2*np.pi,i2):
            return True
        if r>r_i(theta,i1) and r<r_i(theta,i2):
            return True
        if r>r_i(theta+2*np.pi,i1) and r<r_i(theta+2*np.pi,i2):
            return True
        else:
            return False
    if i1==7:
        if r>r_i(theta-2*np.pi,i1) and r<r_i(theta,i2):
            return True
        if r>r_i(theta,i1) and r<r_i(theta+2*np.pi,i2):
            return True
        if r>r_i(theta+2*np.pi,i1) and r<r_i(theta+4*np.pi,i2):
            return True
        else:
            return False
    
def L(Z,H,W):                                                 #disk-halo transition function
    return (1.0+np.e**(-2.0*(abs(Z)-H)/W))**-1.0
    
def regiondebug(r,theta): #outputs region number, 0 to 7.
    regionsfound=0
    region_numbers=False
    
    for n in range(8):
        if regioncheck(r,theta,n,(n+1)%8):
            regionsfound+=1
            region_numbers=(n,(n+1)%8)
            
    return regionsfound, region_numbers
        

#=======X-FIELD PARAMATERS & FUNCTIONS=======
B_X=4.6
Theta0_X=np.radians(49.0)
rc_X=4.8
r_X=2.9

def R_p(R,Z):
    if abs(Z) >= np.tan(Theta0_X)*(R-rc_X):
        return R*rc_X/(rc_X+abs(Z)/np.tan(Theta0_X))

    else:
        return R-abs(Z)/np.tan(Theta0_X)

def b_Xf(R_P):
    return B_X*np.e**(-R_P/r_X)

#def Theta_X(R,Z):
#    if abs(Z)>=np.tan(Theta0_X)*np.sqrt(R)-np.tan(Theta0_X)*rc_X:
#        return np.arctan(abs(Z)/(R-r_p(R,Z)))  #think about atan and maybe problems with quadrants
#    else:
#        return Theta0_X

#striation parameters (the striated field is currently unused)
#gamma= 2.92
#alpha= 2.65 #taken from page 9, top right paragraph
#beta= gamma/alpha-1

#=======OTHER PRELIMINARY DEFINITIONS=======

def mag(V):
    return np.sqrt(V[0]**2+V[1]**2+V[2]**2)
    
def Theta(X,Y): #this theta goes from 0 to 2pi
    if X == 0.0:
        if Y > 0.0:
            return np.pi/2
        if Y < 0.0:
            return 3.0*np.pi/2
        else:
            return 0
    else:
        if Y < 0.0:
            return np.arctan2(Y,X)+2*np.pi
        else:
            return np.arctan2(Y,X)
            
    
def Gam(V):
    return (1.0-(mag(V)/c)**2)**(-0.5)
    
    #check the KE
def KE(VEL):
    if Gam(VEL) != 1:
        return m*c**2*(Gam(VEL)-1)*5.942795e48     #KE, converted to GeV
    else:
        return m*mag(VEL)**2*5.942795e48/2
        
        

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
    
def get_B(pos):
    r=np.sqrt(pos[0]**2.0+pos[1]**2.0)
    theta=Theta(pos[0],pos[1])
    r_p=R_p(r,pos[2])
    b_X=b_Xf(r_p)
    bfield=np.array([0.0, 0.0, 0.0]) #0 for r>1, <20
    
    if mag(pos) > 1.0 and r<20 and r!=0:
        phi_hat=np.array([-pos[1]/r, pos[0]/r, 0.0])
        
#        #halo:
        if pos[2] >= 0.0:
            bfield = np.e**(-abs(pos[2])/z_0)*L(pos[2],h_disk,w_disk)* B_n*(1-L(r,r_n,w_h)) * phi_hat

        else:
            bfield = np.e**(-abs(pos[2])/z_0)*L(pos[2],h_disk,w_disk)* B_s*(1-L(r,r_s,w_h)) * phi_hat
            
        
#        #X-field: 
        if pos[2] != 0: 
            
            bhat_X=(pos[2]**2+(r-r_p)**2)**(-0.5)*np.array([pos[0]*(r-r_p)/r, pos[1]*(r-r_p)/r, pos[2]])
            if pos[2]<0:
                bhat_X*=-1
                
            if abs(r_p) <  rc_X:
                bfield += b_X*(r_p/r)**2*bhat_X

            if abs(r_p) >= rc_X:
                bfield += b_X*(r_p/r)*bhat_X
        else:
            bfield +=b_X*np.array([0,0,1])
                
        #disk:
        if r>= 3.0 and r < 5.0 :
            bfield+=b_ring*(1-L(pos[2],h_disk,w_disk))*phi_hat
            
        if r>=5.0:
            if regioncheck(r,theta,7,0): #region 1
                bfield+=(b1/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([np.sin(theta+i),-np.cos(theta+i),0.0])
                
            elif regioncheck(r,theta,0,1): #region 2
                bfield+=(b2/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([np.sin(theta+i),-np.cos(theta+i),0.0])
                
            elif regioncheck(r,theta,1,2): #region 3
                bfield+=(b3/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([np.sin(theta+i),-np.cos(theta+i),0.0])
                
            elif regioncheck(r,theta,2,3): #region 4
                bfield+=(b4/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([np.sin(theta+i),-np.cos(theta+i),0.0])
                
            elif regioncheck(r,theta,3,4): #region 5
                bfield+=(b5/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([np.sin(theta+i),-np.cos(theta+i),0.0])
                
            elif regioncheck(r,theta,4,5): #region 6
                bfield+=(b6/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([np.sin(theta+i),-np.cos(theta+i),0.0])
                
            elif regioncheck(r,theta,5,6): #region 7
                bfield+=(b7/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([np.sin(theta+i),-np.cos(theta+i),0.0])
                
            elif regioncheck(r,theta,6,7): #region 8
                bfield+=(b8/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([np.sin(theta+i),-np.cos(theta+i),0.0])
                
            elif pos[2]==2:
                print("=====NO REGION FOUND! BUG IN PROGRAM!!!=====")
                print("pos:",pos, end="\n \n ")
                try:
                    global posbug_list
                    posbug_list+=[[x for x in pos]]
                except NameError:
                    posbug_list=[[x for x in pos]]
                    global posbug_list
                    
                    
        bfield=bfield*microGtoT  #convert to Tesla
    return bfield

#FIND THE MAXIMUM B MAGNITUDE ON An XY CROSS SECTION (debug)
#maxB=0.0
#for xxx in np.linspace(-20,20,num=200):
#    for yyy in np.linspace(-20,20,num=200):
#        if xxx**2+yyy**2<25:
#            if mag(get_B([xxx,yyy,.01]))>maxB:
#                maxB=mag(get_B([xxx,yyy,.01]))
#print("MAXB IS",maxB)

def acc_relativistic(pos,vel):   #returns relativistic acceleration. only necessary for large velocity.
    bfield=get_B(pos)
    gam=Gam(vel)
    force = q_b*bfield    #In kg*kpc/s^2
    acc=(gam*m*c**2)**(-0.5)*(c**2*force - np.dot(force,vel)*vel)
    return acc
    
def acc_classical(pos):
    bfield=get_B(pos)
    acc  = bfield*(q_b/m)
    
    return acc
    
#idea: make a relativistic and nonrel paradigm. when the speed is nonrel (difference
#would be less than machine precision) relativistic=False. Every X iterations, check if 
#the paradigm should be changed, where X is the least number of iterations it would take
#to change the velocity by, say, 1% of the speed of light.


tzero=gettimesecs()

g=open(nameofscript+"_full_output"+str(gettimestd())+".txt","a")

#=======PROGRAM BEGINS=======

#for V in [300*(3.24077929e-20)]:
#    for THETA in [-np.pi/2, 0, np.pi/2]:
#        for PHI in [0, np.pi/2, np.pi, 3*np.pi/2]:
tstartsecs=gettimesecs()
tstartstd=gettimestd()

print("=======started at",tstartstd,"with V,THETA,PHI as",V,THETA,PHI,"=======\n")

a=nameofscript+gettimestd()+".txt"
f = open(a, 'a')


vel_0=V*np.array([np.sin(THETA)*np.cos(PHI), np.sin(THETA)*np.sin(PHI), np.cos(THETA)])

pos=np.array([0.0, 0.0, 0.0])
vel=np.array([0.0, 0.0, 0.0])

for n in [0,1,2]:
    pos[n], vel[n]=pos_0[n], vel_0[n]

f.write("\n=====STARTING CONDITIONS FOR RUN "+str(V)+str(THETA)+str(PHI)+" ===== vel_0:\n")
f.write(np.array_str(vel_0)+"\n")
f.write("mag(vel_0):\n")
f.write(str(mag(vel_0))+"\n")
f.write("mag(vel_0)/c\n")
f.write(str(mag(vel_0)/c)+"\n")
f.write("THETA_0\n")
f.write(str(THETA)+"\n")
f.write("PHI_0\n")
f.write(str(PHI)+"\n")
f.write("q_b\n")
f.write(str(q_b)+"\n")
f.write("m\n")
f.write(str(m)+"\n")
f.write("dt\n")
f.write(str(dt)+"\n")
f.write("distance to halt\n")
f.write(str(distance_to_halt)+"\n")

maxvelocity=mag(vel)
distance_tracked=0.0                  #set the distance travelled so far to 0
time=0.0
iterations=0
arclength_in_rho1kpc_region=0.0

#FOR GRAPHING:
xlist, ylist, zlist=[np.array([]) for _ in range(3)]

while mag(pos) < 20.0 and distance_tracked < distance_to_halt:
    iterations+=1
    
    acc=acc_relativistic(pos,vel)
    k1=np.array([-vel*dt,                    -acc_relativistic(pos,           vel)*dt           ])
    k2=np.array([-(vel+k1[1]/2.0)*dt,        -acc_relativistic(pos+k1[0]/2.0, vel+k1[1]/2.0)*dt ])
    k3=np.array([-(vel+k2[1]/2.0)*dt,        -acc_relativistic(pos+k2[0]/2.0, vel+k2[1]/2.0)*dt ])
    k4=np.array([-(vel+k3[1])*dt,            -acc_relativistic(pos+k3[0],     vel+k3[1])*dt     ])
    
    pos_step=(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6.0
    
    #tested: mutability of np array does not mess this up (changing vel, pos does not chg k1)
    vel+=(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6.0
    
    pos+=pos_step
    distance_tracked += mag(pos_step)
    
    #FOR GRAPHING:
    xlist=np.concatenate((xlist,np.array([pos[0]])))
    ylist=np.concatenate((ylist,np.array([pos[1]])))
    zlist=np.concatenate((zlist,np.array([pos[2]])))
    
    
    
    #here, we adjust the timestep so that: 
    #1. the velocity changes by .001% every iteration (hence the 0.00001)
    #2. the change in distance is never greater than 0.01kpc in one iteration (hence the 0.01). equation comes from v^2=v_0^2 + 2a(dx)
#    dt=min([1e50, abs((1/mag(acc))*(mag(vel)+np.sqrt(mag(vel)**2+2*mag(acc)*0.01))) ]) #.000001*mag(vel)/mag(acc)
#    print("dt is",dt)
#    print("dt choices were:",.000001*mag(vel)/mag(acc), abs((1/mag(acc))*(mag(vel)+np.sqrt(mag(vel)**2+2*mag(acc)*0.01))))
#    print()
#    print("pos mag is ",mag(pos))
#    print("vel mag/c is",mag(vel)/c)
#    print("acc mag/c is",mag(acc)/c)
#    print()
    
    #EULER METHOD
#     pos += -vel*dt - 0.5*acc*dt**2
#     
#     vel += -acc*dt   #this is the time-forwards velocity
    
    if mag(vel)>c:
        g.open("SPEED ERROR "+a,"a")
        g.write("beta was "+str(mag(vel)/c))
        print("SPEED ERROR. beta was "+str(mag(vel)/c))
        raise SystemExit()
    
    if iterations%20000==1:
        print("============done with iteration no.",iterations,"============")
        print("pos: ", pos,"\n")
        
        print("vel/c: ",vel/c)
        print("|vel/c|: ",mag(vel)/c,"\n")
        
        print("acc: ", acc)
        print("|acc|: ",mag(acc),"\n")
        
        print("bfield: ",get_B(pos))
        print("|bfield|: ",mag(get_B(pos)))
        print()
        print("displacement: ", mag(pos-pos_0),"kpc")
        print("arc length traversed: ",distance_tracked,"kpc")
        print("simulated time: ",-iterations*dt,"seconds.")
        print()
        print("runtime so far: ",gettimesecs()-tstartsecs,"real seconds")
        print("\n \n")
        
    if mag(pos)<1:
        arclength_in_rho1kpc_region+=pos_step
    
distance_from_start=np.sqrt( (pos[0]-pos_0[0])**2 +(pos[1]-pos_0[1])**2 +(pos[2]-pos_0[2])**2)
theta_f=np.arccos(pos[2]/mag(pos))
phi_f = Theta(pos[0],pos[1])

if mag(pos)<20:
    exitstatus=False

f.write("===FINAL CONDITIONS=== vel:\n")
f.write(np.array_str(vel)+"\n")
f.write("mag(vel)\n")
f.write(str(mag(vel))+"\n")
f.write("mag(vel)/c\n")
f.write(str(mag(vel)/c)+"\n")
f.write("pos\n")
f.write(np.array_str(pos)+"\n")
f.write("mag(pos)\n")
f.write(str(mag(pos))+"\n")
f.write("theta_f\n")
f.write(str(theta_f)+"\n")
f.write("phi_f\n")
f.write(str(phi_f)+"\n")
f.write("distance tracked\n")
f.write(str(distance_tracked)+"\n")
f.write("distance from start\n")
f.write(str(distance_from_start)+"\n")
f.write("Kinetic Energy NONRELATIVISTIC\n")
f.write(str(0.5*m*mag(vel)**2)+"\n")
f.write("maxvelocity\n")
f.write(str(maxvelocity)+"\n")
f.write("maxvelocity/c\n")
f.write(str(maxvelocity/c)+"\n")
f.write("arclength_in_rho1kpc_region\n")
f.write(str(arclength_in_rho1kpc_region)+"\n")
f.write("time\n")
f.write(str(-iterations*dt)+"\n")
f.write("iterations\n")
f.write(str(iterations)+"\n")
f.write("real runtime\n")
f.write(str(gettimesecs()-tstartsecs)+"\n")
f.write("final acc\n")
f.write(np.array_str(acc)+"\n\n\n")
f.write("exit status\n")
f.write(str(exitstatus))
f.close()
g.write("theta, phi, velx, vely, velz for run: "+str(V)+"   "+str(THETA)+"   "+str(PHI)+":\n \n")
g.write(str(theta_f)+"\n")
g.write(str(phi_f)+"\n")
g.write(str(vel[0])+"\n")
g.write(str(vel[1])+"\n")
g.write(str(vel[2])+"\n \n")

g.close()

print("done with run (V,THETA,PHI):  "+str(V)+"   "+str(THETA)+"   "+str(PHI))
print("final speed over c:",mag(vel)/c)
print("iterations:",iterations)
print("done with all runs")
print("finished at",str(gettimestd()))
print("Total time running:",gettimesecs()-tzero,"seconds.")


diskplotbool  =False
trajectorybool=True
vecplotbool   =True
xycrossbool   =False
xzcrossbool   =False

execfile("grapher.py")
