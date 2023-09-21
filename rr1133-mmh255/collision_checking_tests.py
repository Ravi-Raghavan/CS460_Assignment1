import numpy as np 
import matplotlib.pyplot as plt
from collision_checking import *

polygons = np.load("rr1133-mmh255/collision_checking_polygons.npy", allow_pickle= True)
fig, ax = plt.subplots(dpi = 100)
ax.set_aspect("equal")

for index1 in range(len(polygons)):
    polygon = polygons[index1]
    collision_polygon = False
    for index2 in range(len(polygons)):
        if (index1 == index2):
            continue
        if (collides(polygon, polygons[index2])):
            collision_polygon = True
            break
    
    if collision_polygon:
        ax.fill([vertex[0] for vertex in polygon], [vertex[1] for vertex in polygon], alpha=.25, fc='gray', ec='black')
    else:
        ax.fill([vertex[0] for vertex in polygon], [vertex[1] for vertex in polygon], alpha=.25, fc='white', ec='black')


plt.show()