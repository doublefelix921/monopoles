import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

x=2
v=0
dt=0.1
steps=0
bsteps=0

def a(X):
    return X/10+1

while steps < 10000:
    x=x+v*dt+0.5*a(x)*dt**2
    steps+=1
    
print(x)

while bsteps<10000:
    x=x-v*dt-0.5*a(x)*dt**2
    bsteps+=1
    
print(x)