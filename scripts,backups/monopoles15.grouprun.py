#Tracking Magnetic Monopoles Through the Galaxy:
#input starting conditions into the program and it will track the magnetic monopole's path through the bfield.

#Bfield given by "A New Model of the Galatic Magnetic Field" Jansson, Farrar (2012).
#Mass bounds from "Signatures for a Cosmic Flux of Magnetic Monopoles" Wick (2002) are 40TeV <~ M <~ 10^8TeV.
#Wick(2002) has more information on likely mass values.

#units used are all SI units except: distance is in kpc, and angles are in degrees, magnetic field strength is in microgauss.

#variable names generally match those in the Bfield paper.

import math
import numpy as np
import sys
import random
import datetime

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
    
tzero=gettimesecs()

g=open(nameofscript+"_full_output"+str(gettimestd())+".txt","a")

#=======PROGRAM BEGINS=======

for V in [300*(3.24077929e-20)]:
    for THETA in [-math.pi/2, 0, math.pi/2]:
        for PHI in [0, math.pi/2, math.pi, 3*math.pi/2]:
            tstartsecs=gettimesecs()
            tstartstd=gettimestd()
            
            print("=======started at",tstartstd,"with V,THETA,PHI as",V,THETA,PHI,"=======\n")
            
            a=nameofscript+gettimestd()+".txt"
            f = open(a, 'a')
            
            vel_0=V*np.array([math.sin(THETA)*math.cos(PHI), math.sin(THETA)*math.sin(PHI), math.cos(THETA)])
            pos_0=np.array([posx,posy,posz])
            
            pos=np.array([0.0, 0.0, 0.0])
            vel=np.array([0.0, 0.0, 0.0])
            for n in [0,1,2]:
                pos[n]=pos_0[n]
                vel[n]=vel_0[n]
            
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
            
            while mag(pos) < 20.0 and distance_tracked < distance_to_halt and mag(vel) < c:
                r=math.sqrt(pos[0]**2.0+pos[1]**2.0)
                theta=Theta(pos[0],pos[1])
                r_p=R_p(r,pos[2])
                b_X=b_Xf(r_p)
                gam=Gam(vel)
                
                
                bfield=np.array([0.0, 0.0, 0.0]) #0 for r>1, <20
                if mag(pos) > 1.0:
                    #halo:
                    if pos[2] >= 0.0:
                        bfield = math.e**(-abs(pos[2])/z_0)*L(pos[2],h_disk,w_disk)* B_n*(1-L(r,r_n,w_h)) * np.array([-1*pos[1]/r, pos[0]/r, 0.0])
                    else:
                        bfield = math.e**(-abs(pos[2])/z_0)*L(pos[2],h_disk,w_disk)* B_s*(1-L(r,r_s,w_h)) * np.array([-1*pos[1]/r, pos[0]/r, 0.0])
                        
                    #X-field: 
                    if pos[2] != 0: 
                        bhat_X=(pos[2]**2+(r-r_p)**2)**(-0.5)*np.array([pos[0]*(r-r_p)/r, pos[1]*(r-r_p)/r, pos[2]])
                        if pos[2] < 0.0:
                            bhat_X[2]=-1.0*bhat_X[2]
                        if abs(r_p) <  rc_X:
                            bfield += b_X*(r_p/r)**2*bhat_X
                        if abs(r_p) >= rc_X:
                            bfield += b_X*(r_p/r)*bhat_X
                            
                    #disk:
                    if r>= 3.0 and r < 5.0 :
                        bfield+=b_ring*(1-L(pos[2],h_disk,w_disk))
                        
                    elif r>=5.0 and r<=20.0:
                        if   r >= r_i(theta,7) and r < r_i(theta,0): #region 1
                            bfield+=(b1/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
                            
                        elif r >= r_i(theta,0) and r < r_i(theta,1): #region 2
                            bfield+=(b2/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
                            
                        elif r >= r_i(theta,1) and r < r_i(theta,2): #region 3
                            bfield+=(b3/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
                            
                        elif r >= r_i(theta,2) and r < r_i(theta,3): #region 4
                            bfield+=(b4/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
                            
                        elif r >= r_i(theta,3) and r < r_i(theta,4): #region 5
                            bfield+=(b5/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
                            
                        elif r >= r_i(theta,4) and r < r_i(theta,5): #region 6
                            bfield+=(b6/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
                            
                        elif r >= r_i(theta,5) and r < r_i(theta,6): #region 7
                            bfield+=(b7/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
                            
                        elif r >= r_i(theta,6) and r < r_i(theta,7): #region 8
                            bfield+=(b8/r)*(1.0-L(pos[2],h_disk,w_disk))*np.array([math.sin(theta+i),-math.cos(theta+i),0.0])
                            
                #NONREL
                #acc  = bfield*(q_b/m)
                #pos -= vel*dt + 0.5*acc*dt**2
                #distance_tracked += math.sqrt((vel[0]*dt+0.5*acc[0]*dt**2)**2+(vel[1]*dt+0.5*acc[1]*dt**2)**2+(vel[2]*dt+0.5*acc[2]*dt**2)**2)
                
                #vel += acc*dt #this is the time-forwards velocity
                
                #RELATIVISTIC
                force = q_b*bfield #In kg*kpc/s^2
                
                acc=(gam*m*c**2)**(-0.5)*(c**2*force - np.dot(force,vel)*vel)
                
                #here, we adjust the timestep so that: 
                #1. the velocity changes by .001% every iteration (hence the 0.00001)
                #2. the change in distance is never greater than 0.01kpc in one iteration (hence the 0.01). equation comes from v^2=v_0^2 + 2a(dx)
                dt=min([1e50, abs((1/mag(acc))*(mag(vel)+math.sqrt(mag(vel)**2+2*mag(acc)*0.01))) ]) #.000001*mag(vel)/mag(acc)
                print("dt is",dt)
                print("dt choices were:",.000001*mag(vel)/mag(acc), abs((1/mag(acc))*(mag(vel)+math.sqrt(mag(vel)**2+2*mag(acc)*0.01))))
                print()
                print("pos mag is ",mag(pos))
                print("vel mag/c is",mag(vel)/c)
                print("acc mag/c is",mag(acc)/c)
                print()
                
                pos += -vel*dt - 0.5*acc*dt**2
                
                distance_tracked += math.sqrt((vel[0]*dt+0.5*acc[0]*dt**2)**2+(vel[1]*dt+0.5*acc[1]*dt**2)**2+(vel[2]*dt+0.5*acc[2]*dt**2)**2)
                
                vel += acc*dt
                
                time-=dt
                
                if mag(vel)>c:
                    g.open("SPEED ERROR "+a,"a")
                    g.write("spd was",str(mag(vel)/c))
                
#                 if random.randint(1,10000)==1:
#                     print("================iteration no.",-time/dt,"================")
#                     print("Position (kpc) is", pos)
#                     print("...with a magnitude of",mag(pos))
#                     print("Velocity (kpc/s) is",vel)
#                     print("...with a magnitude of",mag(vel)/c,"times c")
#                     print("Acceleration(kpc/s/s) is", acc)
#                     print("...with a magnitude of",mag(acc))
#                     print("bfield is",bfield)
#                     print("with a magnitude of",mag(bfield))
#                     print()
#                     print("Distance from initial position is", math.sqrt( (pos[0]-pos_0[0])**2 +(pos[1]-pos_0[1])**2 +(pos[2]-pos_0[2])**2),"kpc")
#                     print("Distance traveled so far is",distance_tracked,"kpc")
#                     print("Current time is ",time,"seconds.")
#                     print()
#                     print("runtime so far is",gettimesecs()-tstartsecs,"real seconds")
#                     print()
                
            distance_from_start=math.sqrt( (pos[0]-pos_0[0])**2 +(pos[1]-pos_0[1])**2 +(pos[2]-pos_0[2])**2)
            theta_f=math.acos(pos[2]/mag(pos))
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
            f.write("time\n")
            f.write(str(time)+"\n")
            f.write("iterations\n")
            f.write(str(-time/dt)+"\n")
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
            print("done with run (V,THETA,PHI):  "+str(V)+"   "+str(THETA)+"   "+str(PHI))
            print("final speed over c:",mag(vel)/c)
            print("iterations:",-time/dt)

print("done with all runs")
print("finished at",str(gettimestd()))
print("Total time running:",gettimesecs()-tzero,"seconds.")

g.close()
