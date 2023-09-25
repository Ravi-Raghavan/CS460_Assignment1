import numpy as np 
import matplotlib.pyplot as plt
from collision_checking import *

#Testing this on multiple polygon scenes
polygons = np.load("rr1133-mmh255/polygons_scene8.npy", allow_pickle= True)

#Pre-process polygon list to numpy array, just-in-case
for index, polygon in enumerate(polygons):
    polygon = [list(vertex) for vertex in polygon]
    polygons[index] = np.array(polygon)

fig, ax = plt.subplots(dpi = 100)
ax.set_aspect("equal")

for index1 in range(len(polygons)):
    polygon = polygons[index1]
    collision_polygon = False
    for index2 in range(len(polygons)):
        if (index1 == index2):
            continue
        if (collides_optimized(polygon, polygons[index2])):
            collision_polygon = True
            break
    
    if collision_polygon:
        ax.fill([vertex[0] for vertex in polygon], [vertex[1] for vertex in polygon], alpha=.25, fc='gray', ec='black')
    else:
        ax.fill([vertex[0] for vertex in polygon], [vertex[1] for vertex in polygon], alpha=.25, fc='white', ec='black')


plt.show()