import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import random
import math

# receive input about scene -> I need to change this to read from a file lol
P = int(input('number of polygons: '))
nmin = int(input('minimum vertices: '))
nmax = int(input('maximum vertices: '))
rmin = float(input('minimum radius: '))
rmax = float(input('maximum radius: '))

# store polygons
polygons = []

# create polygon P times
for _ in range(P):
    # store vertices
    polygon = []

    # get center for polygon
    x_center = random.uniform(0,2)
    y_center = random.uniform(0,2)

    # get N vertices in range [nmin,nmax] for the centers
    N = random.randint(nmin,nmax)
    for _ in range(N):
        # get angle a using range [0,360] & convert to radians
        a_deg = random.uniform(0,360)
        a_rad = math.radians(a_deg)

        # place vertex on a circle (of radius rmin) around the center 
        x_vertex = x_center + rmin*math.cos(a_rad)
        y_vertex = y_center + rmin*math.sin(a_rad)

        # get r using range [0,rmax-rmin] and extrude the vertex
        r = random.uniform(0, rmax-rmin)
        x = x_vertex + r*math.cos(a_rad)
        y = y_vertex + r*math.sin(a_rad)

        # add vertex to polygon
        polygon.append((x,y))
    polygons.append(np.array(polygon))

# check if polygons are convex using ConvexHull()
convexPolygons = []
for i in polygons:
   hull = ConvexHull(i)
   #extract each vertex from hull
   vertices = [i[j] for j in hull.vertices]   
   convexPolygons.append(np.array(vertices))

# plot polygons
fig, ax = plt.subplots(figsize=(7,7), dpi=100)
for i in convexPolygons:
    # save polygon -> figuring out how to save the polygons
    #np.save('polygon_scene.npy', i)
    x_coordinates, y_coordinates = zip(*i)
    #close polygon
    x_coordinates = np.append(x_coordinates, x_coordinates[0])
    y_coordinates = np.append(y_coordinates, y_coordinates[0])

    ax.plot(x_coordinates, y_coordinates)

ax.set(xlim = (0,2), ylim =(0,2))
plt.show()