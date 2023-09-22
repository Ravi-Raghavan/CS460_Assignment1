import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import random
import math

# get user input 
P = int(input('number of polygons: '))
nmin = int(input('minimum vertices: '))
nmax = int(input('maximum vertices: '))
rmin = float(input('minimum radius: '))
rmax = float(input('maximum radius: '))

# construct polygonal map
def constructScene(P):
    polygons = []
    convexPolygons = []

    # construct P polygons
    for _ in range(P):
        polygons.append(constructPolygon(nmin, nmax, rmin, rmax))

    # check each polygon is convex 
    for i in polygons:
        makeConvex(i, convexPolygons)
    
    # save polygons to .npy file
    polygons_save = np.array(convexPolygons,dtype = object)
    np.save('polygons_scene.npy', polygons_save)

    # plot the polygons
    visualize(convexPolygons)

# construct polygon using ranges for vertices and radii
def constructPolygon(nmin, nmax, rmin, rmax):
    polygon =[]
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
    return polygon

# makes a polygon convex using ConvexHull
def makeConvex(polygon, convexPolygons):
    # get indexes of points that make polygon convex
    hull = ConvexHull(polygon)
    # extract each vertex from hull 
    vertices = [polygon[j] for j in hull.vertices]   
    convexPolygons.append(vertices)

# plots the convexPolygon 
def visualize(convexPolygons):
    fig, ax = plt.subplots(figsize=(7,7), dpi=100)
    for i in convexPolygons:
         # store x and y coordinates in seperate arrays
         x_coordinates, y_coordinates = zip(*i)
         # closes polygon
         x_coordinates = np.append(x_coordinates, x_coordinates[0])
         y_coordinates = np.append(y_coordinates, y_coordinates[0])
         # plot points
         ax.plot(x_coordinates, y_coordinates)
    ax.set(xlim = (0,2), ylim =(0,2))
    ax.set_aspect('equal')
    plt.show() 

constructScene(P)