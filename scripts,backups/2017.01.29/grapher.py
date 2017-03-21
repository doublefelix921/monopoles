"""
This program contains incomplete code, and only works when addended to the main monopoles.py file. It creates plots of:
1. The disk regions                           (if diskplotbool  =True)
2. The monopole's trajectory                  (if trajectorybool=True)
3. A 3D vector bfield plot                    (if vecplotbool   =True)
4. An XY cross section of the magnitude of B  (if xycrossbool   =True)
5. An XZ cross section of the magnitude of B  (if xzcrossbool   =True)
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#meta parameters
# if circle3 is True, a radius-5 circle will be plotted centered on the
# galaxy, for reference. This is for the vector plot only.
circle1 = False
circle3 = False
circle5 = False
circle20 = False

#=========PLOTTING THE REGIONS========

if diskplotbool==True:
    print("plotting fig1: disk regions...\n")
    plt.figure(1)
    polar_ax=plt.subplot(111, polar=True)
    polar_ax.set_rlim((0, 20))
    
    #list of regiondata, each region has an rlist and a thetalist
    plotdata=[ [[], []] for n in range(8)]
    
    for rr in np.linspace(5, 20, num=50):
        for ttheta in np.linspace(0, 2*np.pi, num=500):
            howmanyregions=regiondebug(rr, ttheta)[0]
            
            if howmanyregions==1:
                regionnumber=(regiondebug(rr, ttheta)[1][1])%8+1 
                #the region i'm in is the biggest index of the surrounding curves + 1. 
                
                plotdata[regionnumber-1][0]+=[ttheta]
                plotdata[regionnumber-1][1]+=[rr]
    
    colorlist=["bo","go","ro","co","mo","yo","ko","wo"]
    
    for k in range(8):
        polar_ax.plot(plotdata[k][0], plotdata[k][1], colorlist[k], markersize=3)
    
    plt.title("Map of disk field regions")
    
    plt.text(np.pi, 32, 
        "\nblue: r1\ngreen: r2\nred: r3\ncyan: r4\nmagenta:" + 
        "r5\nyellow: r6\nblack: r7\ngray: r8 ")

if trajectorybool==True:
    #=====PLOTTING THE TRAJECTORY=====
    print("plotting fig2: trajectory...\n")
    trajectory_fig=plt.figure(2)
    trajectory_ax=trajectory_fig.gca(projection="3d")
    trajectory_ax.plot(xlist, ylist, zlist, "r")
    
    trajectory_ax.set_title("Trajectory of Monopole")
    
    trajectory_ax.text2D(0.00, 0.00, "x0="+str(r0)+"\nv0="+str(v0), 
        transform=trajectory_ax.transAxes)
    
    trajectory_ax.set_xlabel('X')   #label each positional axis
    trajectory_ax.set_ylabel('Y')
    trajectory_ax.set_zlabel('Z')
    
    trajectory_ax.set_xlim3d(-20, 20)  #set view ranges for the x, y, and z axes  
    trajectory_ax.set_ylim3d(-20, 20)
    trajectory_ax.set_zlim3d(-20, 20)  

if vecplotbool==True:
    #=====PLOTTING THE BFIELD VECTOR FIELD=====
    print("plotting fig3: bfield...\n")
    vector_fig=plt.figure(3)
    vector_ax=vector_fig.gca(projection="3d")
    vector_ax.set_title("Milky Way Bfield")
    xmin, xmax=--2, 2 
    ymin, ymax=-20, 20
    zmin, zmax=-20, 20
    
    arrowlength=1.3
    
    plotstep=3  #what is the spacing between adjacent vectors along axes? 
                #(defines the density of arrowsh)
    x, y, z   = np.meshgrid(np.arange(xmin, xmax+1, plotstep),
                          np.arange(ymin, ymax+1, plotstep),
                          np.arange(zmin, zmax+1, plotstep))
    #for some reason, the above command creates x, y, z
    #matricies with lists of coordinates (y,x,z). I don't
    #know why it would order the coordinates like that.
    #but it's the way that it works in the quiver example 
    #too, and it seems to graph correctly based on their ranges.
    
    howmany_x, howmany_y, howmany_z = (xmax-xmin)//plotstep+1, (ymax-ymin)//plotstep+1, (zmax-zmin)//plotstep+1
    
    
    u=np.zeros(shape=(howmany_y, howmany_x, howmany_z)) #shape is (number of x values being plotted,  num y values, num z values)
    v=np.zeros(shape=(howmany_y, howmany_x, howmany_z))
    w=np.zeros(shape=(howmany_y, howmany_x, howmany_z))
    
    
    def f(VEC):
        return np.array([VEC[0], 0, 0])
        
        
    for I in range(howmany_y):
        for J in range(howmany_x):   
            for K in range(howmany_z):
                position=np.array([x[I][J][K], y[I][J][K], z[I][J][K]])
    #            print(position)
    #            print("")
    #            print("i,j,k are", I,J,K)
    #            print("position is:", position)
    #            print("r is:", mag([position[0],position[1],0]))
                u[I][J][K]=get_B(position)[0]*100  #these are multplied by 100
                v[I][J][K]=get_B(position)[1]*100  #because they otherwise aren't
                w[I][J][K]=get_B(position)[2]*100  #plotted for some reason.
                
    vector_ax.quiver(x, y, z, u, v, w, length=arrowlength)
    
    #PLOTTING CIRCLES OF R=1, 3, 5, 20
    x_rad1_list, y_rad1_list, x_rad3_list, y_rad3_list, x_rad5_list, y_rad5_list, x_rad20_list, y_rad20_list, =[ [] for _ in range(8) ]
    
    for angle in np.linspace(0.0, 2*np.pi, num=int(2*np.pi//0.1)):
        x_rad1_list+=[1*np.cos(angle)]
        y_rad1_list+=[1*np.sin(angle)]
        x_rad3_list+=[3*np.cos(angle)]
        y_rad3_list+=[3*np.sin(angle)]
        x_rad5_list+=[5*np.cos(angle)]
        y_rad5_list+=[5*np.sin(angle)]
        x_rad20_list+=[20*np.cos(angle)]
        y_rad20_list+=[20*np.sin(angle)]
    
    z_circle=[]
    for n in range(len(x_rad1_list)):
        z_circle+=[0]
        
    if circle1==True:
        vector_ax.plot(x_rad1_list, y_rad1_list, z_circle, "b")
    if circle3==True:
        vector_ax.plot(x_rad3_list, y_rad3_list, z_circle, "r")
    if circle5==True:
        vector_ax.plot(x_rad5_list, y_rad5_list, z_circle, "g")
    if circle20==True:
        vector_ax.plot(x_rad20_list, y_rad20_list, z_circle, "m")
    
    vector_ax.set_xlabel('X')   #label each positional axis
    vector_ax.set_ylabel('Y')
    vector_ax.set_zlabel('Z')
    
    # set view ranges for the x, y, and z axes  
    vector_ax.set_xlim3d(xmin-arrowlength, xmax+arrowlength)  
    vector_ax.set_ylim3d(ymin-arrowlength, ymax+arrowlength)
    vector_ax.set_zlim3d(zmin-arrowlength, zmax+arrowlength)     

#=====BUGTESTING PLOT: XY CROSS SECTION WITH MAGNITUDES =====

if xycrossbool == True:
    print("plotting fig4: xy cross section...\n")
    xy_fig=plt.figure(4)
    xy_ax=plt.subplot(111)
    
    zcrossvalue=0.1
    
    plt.title("XY Cross Section, Z="+str(zcrossvalue)+"kpc")
    
    xylistx, xylisty, xylistBmag = [], [], []
    
    for xx in np.linspace(-21, 21, num=801):
        for yy in np.linspace(-21, 21, num=801):
            xylistBmag+=[mag(get_B([xx, yy, zcrossvalue]))]
            xylistx+=[xx]
            xylisty+=[yy]
            
    xybiggestB=max(xylistBmag)
    xysmallestB=min([_ for _ in xylistBmag if _!=0])
    for n in range(len(xylistBmag)):
        
        if xylistBmag[n] != 0:
            #here, I scale the B magnitudes so that the lowest one has size  ~10^-2 and the greatest one is 1. 
            #Then I flip them about 0.5, according to matplotlib's convention, so that the strongest is black.
            #I also make them a string to be interpretable in matplotlib's color thing.
            xylistBmag[n]=str(1-np.power(xylistBmag[n]/xybiggestB, 2.0/17.0))
        else:
            #zero spots are colored green to distinguish from near-white nonzero spots.
            xylistBmag[n]="g"
    
    
    
    xy_ax.scatter(xylistx, xylisty, color=xylistBmag, marker='o')
    
    xy_ax.set_xbound(-20, 20)
    xy_ax.set_ybound(-20, 20)
    
    xy_ax.set_xlabel("X")
    xy_ax.set_ylabel("Y")
    
    #plt.text(-26,-24,"White: "+str(xysmallestB)+"T. Blk: "+str(xybiggestB)+"T, Green=0T.")
    plt.text(-26, -24, ("White: %.2e" % xysmallestB) + 
        ("T  Black: %.2e" % xybiggestB) + "T, Green=0T.")
if xzcrossbool==True:

    #=====BUGTESTING PLOT: XZ CROSS SECTION WITH MAGNITUDES  =====
    print("plotting fig5: xz cross section...\n")
    xz_fig=plt.figure(5)
    xz_ax=plt.subplot(111)
    
    ycrossvalue=0.1
    
    plt.title("XZ Cross Section, Y="+str(ycrossvalue)+"kpc")
    
    xzlistx, xzlistz, listBmag = [], [], []
    
    xdomain = [-21, 21]
    zdomain = [-51, 51]

    for xx in np.linspace(xdomain[0], xdomain[1], num=401):
        for zz in np.linspace(zdomain[0], zdomain[1], num=401):
            listBmag+=[mag(get_B([xx, ycrossvalue, zz]))]
            xzlistx+=[xx]
            xzlistz+=[zz]
            
    biggestB=max(listBmag)
    smallestB=min([_ for _ in listBmag if _!=0])

    for n in range(len(listBmag)):
        
        if listBmag[n] != 0:
            #here, I scale the B magnitudes so that the lowest one has exponent 10^-2 and the greatest one is 1. 
            #Then I flip them about 0.5, according to matplotlib's convention, so that the strongest is black.
            #I also make them a string to be interpretable in matplotlib's color thing.
            listBmag[n]=str(1-np.power(listBmag[n]/biggestB, 2.0/17.0))
        else:
            #zero spots are colored green to distinguish from near-white nonzero spots.
            listBmag[n]="g"
    
    
    xz_ax.scatter(xzlistx, xzlistz, color=listBmag, marker='o')
    
    xz_ax.set_xbound(xdomain[0], xdomain[1])
    xz_ax.set_ybound(zdomain[0], zdomain[1])
    
    xz_ax.set_xlabel("X")
    xz_ax.set_ylabel("Z")
    
    plt.text(-26, -24, ("White: %.2e" % smallestB) + 
        ("T  Black: %.2e" % biggestB) + "T, Green=0T.")

#=====DONE WITH PLOTS=====
plt.show()
