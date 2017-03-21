import math
import numpy as np

x=-4.2
y=-16.60942886
z=-7.06208855

theta= np.arctan(y/x)
phi=np.arccos(z/np.sqrt(x**2+y**2+z**2))
print(theta, phi)