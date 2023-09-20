import numpy as np 
import matplotlib.pyplot as plt
from collision_checking import *

polygons = np.load("rr1133-mmh255/collision_checking_polygons.npy", allow_pickle= True)

# for index in range(len(polygons)):
#     print(f"INDEX: {index}")
#     print(polygons[index])
#     print("-----------------")

fig, ax = plt.subplots()

for index1 in range(len(polygons)):
    polygon = polygons[index1]
    collision_polygon = False
    for index2 in range(len(polygons)):
        if (index1 == index2):
            continue
        if (collides(polygon, polygons[index2])):
            collision_polygon = True
            if (index1 == 18):
                print(f"WTF: {index2}")
            break
    
    if collision_polygon:
        ax.fill([vertex[0] for vertex in polygon], [vertex[1] for vertex in polygon], alpha=.25, fc='gray', ec='black')
    else:
        ax.fill([vertex[0] for vertex in polygon], [vertex[1] for vertex in polygon], alpha=.25, fc='white', ec='black')


# ax.fill([vertex[0] for vertex in polygons[18]], [vertex[1] for vertex in polygons[18]], alpha=.25, fc='gray', ec='black')
# ax.fill([vertex[0] for vertex in polygons[0]], [vertex[1] for vertex in polygons[0]], alpha=.25, fc='gray', ec='black')

plt.show()

# print(point_contained(polygons[2], np.array([0.17007955, 1.89843194])))
# print(collides(polygons[18], polygons[0]))