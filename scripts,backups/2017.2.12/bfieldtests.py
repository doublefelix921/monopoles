from __future__ import division, absolute_import, print_function, unicode_literals

import numpy as np
import random as r
import usefulfunctions as f
import getB as gb

#this program is called by monopoles.py if it contains bfieldtests == true. 
mag = f.mag

#get avg B (and mean deviation B) over 10,000 random points.
Blist = []
numPoints = 1000000
print("Sampling %d random Bfield points to find min, max, avg, and mean absolute dev..." % numPoints)
for n in range(numPoints):
    probePhi = r.uniform(0, 2*np.pi)
    probeTheta = r.uniform(-np.pi, np.pi)
    probeR = r.uniform(0.0, 20.0)
    probeB = mag(gb.get_B([probeR*np.sin(probeTheta)*np.cos(probePhi), probeR*np.sin(probeTheta)*np.sin(probePhi), probeR*np.cos(probeTheta)]))
    Blist += [probeB]

avg = sum(Blist)/len(Blist)
avgdev = sum([abs(B - avg) for B in Blist])/len(Blist)

print("\nTested 1,000,000 random |B| over the r=20kpc sphere.")
print("avg |B| was", avg, "T +-",avgdev,"(avg abs deviation)")
print("avg |B| was %.4gT +-%.4gT (avg abs deviation)" % (avg, avgdev))

Blist = [B for B in Blist if B!=0]
print("nonzero |B| ranged from%.4gT to %.4gT" % 
    (min(Blist), max(Blist)))