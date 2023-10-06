import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from scipy.spatial import ConvexHull
import random

# construct polygonal map 
def constructScene(P, nmin, nmax, rmin, rmax, name):
    polygons = []
    convexPolygons = []
    
    # construct P polygons
    for _ in range(P):
        contained = True
        # while the polygon is NOT within (0,2), create again
        while contained:
            inBounds = True
            # creates polygon and checks if vertices are within (0,2)
            polygon = constructPolygon(nmin, nmax, rmin, rmax)
            for vertex in polygon: 
                if(vertex[0] >= 2 or vertex[1] >= 2 or vertex[0] <= 0 or vertex[1] <= 0):
                    inBounds = False
            # if polygon is within boundaries, add to polygons
            if inBounds:
                polygons.append(polygon)
                contained = False;

    # check each polygon is convex 
    for polygon in polygons:
        makeConvex(polygon, convexPolygons)

    # save polygons to .npy file
    polygons_save = np.array(convexPolygons, dtype=object)
    np.save(name, polygons_save)

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
        a_rad = np.deg2rad(a_deg)

        # place vertex on a circle (of radius rmin) around the center 
        x_vertex = x_center + rmin*np.cos(a_rad)
        y_vertex = y_center + rmin*np.sin(a_rad)

        # get r using range [0,rmax-rmin] and extrude the vertex
        r = random.uniform(0, rmax-rmin)
        x = x_vertex + r*np.cos(a_rad)
        y = y_vertex + r*np.sin(a_rad)

        # add vertex to polygon
        polygon.append((x,y))
    return polygon

# makes a polygon convex using ConvexHull
def makeConvex(polygon, convexPolygons):
    # get indexes of points that make polygon convex
    hull = ConvexHull(polygon)
    # extract each vertex from hull 
    vertices = [polygon[j] for j in hull.vertices]  
    vertices.append(vertices[0]) 
    vert = np.array(vertices)
    convexPolygons.append(vert)

# plots the convexPolygon 
def visualize(convexPolygons):
    fig, ax = plt.subplots(figsize=(7,7), dpi=100)    
    for polygon in convexPolygons:
        ax.add_patch(Polygon(polygon, closed = True, ec = 'black', facecolor = 'grey', alpha = 0.2))
    ax.set(xlim = (0,2), ylim =(0,2))
    ax.set_aspect('equal')
    plt.show() 

def get_input():
    # get user input 
    P = int(input('number of polygons: '))
    nmin = int(input('minimum vertices: '))
    nmax = int(input('maximum vertices: '))
    rmin = float(input('minimum radius: '))
    rmax = float(input('maximum radius: ')) 
    name = (input('output file name: '))

    constructScene(P,nmin,nmax,rmin,rmax,name)

get_input()