from __future__ import division, print_function
import random
import numpy as np
import usefulfunctions as f
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

try: input = raw_input
except NameError: pass

mag = f.mag

def fibonacci_sphere(samples=1, randomize=True):
    # this function returns many ~evenly distributed points on a unit sphere.
    # it is taken from stackoverflow:
    # http://stackoverflow.com/a/26127012/3347826
    # method explained on http://blog.marmakoide.org/?p=1
    rnd = 1.
    if randomize:
        rnd = random.random() * samples

    # if samples <= 14:
    #     rnd = 0.5 * samples

    points_rectangular = []
    offset = 2./samples
    increment = np.pi * (3. - np.sqrt(5.)); #golden ratio

    for i in range(samples):
        y = ((i * offset) - 1) + (offset / 2);
        r = np.sqrt(1 - pow(y,2))

        phi = ((i + rnd) % samples) * increment

        x = np.cos(phi) * r
        z = np.sin(phi) * r

        points_rectangular.append([x, y, z])

    return points_rectangular
# samp = int(round(random.random()*50))
# a = fibonacci_sphere(samples=samp, randomize=True)
# b = fibonacci_sphere(samples=samp, randomize=True)
# avgarea = 4*np.pi/len(a)
# avgdist = np.sqrt(avgarea)
# avgdist2 = np.sqrt(avgarea/np.pi)
# print("Estimates for expected distance btw pts:\n %.4g, %.4g" % (avgdist, avgdist2))
# print()

# PLOTTING THE DIRECTIONS
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# xs = [p[0] for p in points_rectangular]
# ys = [p[1] for p in points_rectangular]
# zs = [p[2] for p in points_rectangular]
# ax.scatter(xs, ys, zs)
# ax.set_xlabel('X (kpc)')
# ax.set_ylabel('Y (kpc)')
# ax.set_zlabel('Z (kpc)')
# plt.show()

def metrics(pts):
    pts = np.asarray(pts)
    

    mindist = 2
    maxdist = 0

    highestx = 0
    highesty = 0
    highestz = 0
    nearestNeighborDistances = []
    for n in range(len(pts)):
        p = pts[n]
        if p[0]>highestx: highestx = p[0]
        if p[1]>highesty: highesty = p[1]
        if p[2]>highestz: highestz = p[2]
        nnDist = 2
        for m in range(len(pts)):
            if m == n: continue
            q = pts[m]

            dist = mag(p-q)
            if dist < mindist: mindist = dist
            if dist < nnDist: nnDist = dist
            if dist > maxdist: maxdist = dist
            # if dist < avgdist/2:
            #     print("points p and q had a very low avg dist:")
            #     print(p, q)
            #     print(mag(p-q))
            #     print()
        nearestNeighborDistances += [nnDist]
    # avgnnDist = sum(nearestNeighborDistances)/len(nearestNeighborDistances)
    # minnnDist_norm = 
    # print("min, max nearest neighbor distance:\n %.4g, %.4g  *avgdist" % 
    #     (min(nearestNeighborDistances)/avgdist, max(nearestNeighborDistances)/avgdist))
    # # print("highest x,y,z: %.4g, %.4g, %.4g" % (highestx, highesty, highestz))
    # print("Avg nearest neighbor distance:\n %.4g +- %.4g" % (avgnnDist, 
        # np.std(nearestNeighborDistances)))
    # print()
    return nearestNeighborDistances


sampleNums = range(5,101)
avgavgNNDs = [ -1 for n in range(101)]
avgminNNDs = [ -1 for n in range(101)]
avgmaxNNDs = [ -1 for n in range(101)]
avgstdevPerMeans = [ -1 for n in range(101)]
avgNNDs_nonrandom = [ -1 for n in range(101)]
minNNDs_nonrandom = [ -1 for n in range(101)]
maxNNDs_nonrandom = [ -1 for n in range(101)]
stdevPerMeans_nonrandom = [ -1 for n in range(101)]

for sampN in sampleNums:
    # generate fib 20 times randomly
    numRandRuns = 20
    avgNNDs_normed, minNNDs_normed, maxNNDs_normed, stdevPerMeans = [], [], [], []
    avgdist = np.sqrt(4*np.pi/sampN)
    for n in range(numRandRuns):
        NNDs = metrics(fibonacci_sphere(samples = sampN, randomize = True))
        avgNNDs_normed += [np.average(NNDs)/avgdist]
        minNNDs_normed += [min(NNDs)/avgdist]
        maxNNDs_normed += [max(NNDs)/avgdist]
        stdevPerMeans += [np.std(NNDs)/np.average(NNDs)]

    avgavgNNDs[sampN] = np.average(avgNNDs_normed)
    avgminNNDs[sampN] = np.average(minNNDs_normed)
    avgmaxNNDs[sampN] = np.average(maxNNDs_normed)
    avgstdevPerMeans[sampN] = np.average(stdevPerMeans)

    # generate fib once nonrandomly 
    NNDs = metrics(fibonacci_sphere(samples = sampN, randomize = False))
    avgNNDs_nonrandom[sampN] = np.average(NNDs)/avgdist
    minNNDs_nonrandom[sampN] = min(NNDs)/avgdist
    maxNNDs_nonrandom[sampN] = max(NNDs)/avgdist
    stdevPerMeans_nonrandom[sampN] = np.std(NNDs)/np.average(NNDs)



avgmin_minus_avgmax_NNDs = [avgmaxNNDs[n]-avgminNNDs[n] for n in range(len(avgminNNDs))]
min_minus_max_NNDs_nonrandom = [maxNNDs_nonrandom[n]-minNNDs_nonrandom[n] for n in range(len(avgminNNDs))]

# smooth out the graphs
def smoothpoints(points, start):
    for n in range(3):
        for n in range(start, len(points)-1):
            points[n] = 0.5*(points[n]+points[n+1])
    return points

avgavgNNDs = smoothpoints(avgavgNNDs, 5)
avgminNNDs = smoothpoints(avgminNNDs, 5)
avgmaxNNDs = smoothpoints(avgmaxNNDs, 5)
avgstdevPerMeans = smoothpoints(avgstdevPerMeans, 5)
avgNNDs_nonrandom = smoothpoints(avgNNDs_nonrandom, 5)
minNNDs_nonrandom = smoothpoints(minNNDs_nonrandom, 5)
maxNNDs_nonrandom = smoothpoints(maxNNDs_nonrandom, 5)
stdevPerMeans_nonrandom = smoothpoints(stdevPerMeans_nonrandom, 5)


print("plotting average distances graph")
plt.plot(sampleNums, avgavgNNDs[5:101], 'r-')
plt.plot(sampleNums, avgNNDs_nonrandom[5:101], 'b-')
plt.show()
ans1 = input("Press enter when ready for mins graph")
plt.plot(sampleNums, avgminNNDs[5:101], 'r-')
plt.plot(sampleNums, minNNDs_nonrandom[5:101], 'b-')
plt.show()

ans2 = input("Press enter when ready for maxs graph")
plt.plot(sampleNums, avgmaxNNDs[5:101], 'r-')
plt.plot(sampleNums, maxNNDs_nonrandom[5:101], 'b-')
plt.show()

ans3 = input("Press enter when ready for max minus min graph")
plt.plot(sampleNums, avgmin_minus_avgmax_NNDs[5:101], 'r-')
plt.plot(sampleNums, min_minus_max_NNDs_nonrandom[5:101], 'b-')
plt.show()

ans4 = input("Press enter when ready for stdev per mean graph")
plt.plot(sampleNums, avgstdevPerMeans[5:101], 'r-')
plt.plot(sampleNums, stdevPerMeans_nonrandom[5:101], 'b-')
plt.show()
