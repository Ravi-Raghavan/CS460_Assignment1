import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from collision_checking import *


#Print Matplotlib version for sanity
print(f"matplotlib version: {matplotlib.__version__}")

#Initialize figures and axes
f,ax = plt.subplots(dpi = 100)
ax.set_aspect("equal")

#set axis limits for x and y axes
ax.set_xlim(-0.1, 2.2)
ax.set_ylim(-0.1, 2.2)

class MinkowskiPlot:
    def __init__(self, f, ax, rotation_angle):
        #Store figure and axes as instance variables
        self.f = f
        self.ax = ax
        self.rotation_angle = rotation_angle
        
        #Load Polygon Data
        self.polygons = np.load("rr1133-mmh255/2d_rigid_body.npy", allow_pickle= True)

        #Generate Starting Configuration for Rectangle with configuration at origin
        self.rectangle = self.random_rectangle_configuration(self.generate_rectangle())        
             
    #Rotate and Translate the Default Rectangle by a given angle and translation amount
    def random_rectangle_configuration(self, rectangle):
        self.center_point = np.array([0, 0])
        angle = np.deg2rad(self.rotation_angle)
        rotation_matrix = np.array([[np.cos(angle), -1 * np.sin(angle)], [np.sin(angle), np.cos(angle)]])    
        return (rotation_matrix @ rectangle.T).T + self.center_point

    #Generate Rectangle that is centered around the origin
    def generate_rectangle(self):
        w = 0.2
        h = 0.1
        
        top_left = np.array([-1 * w/2, h/2])
        top_right = np.array([w/2, h/2])
        bottom_left = np.array([-1 * w/2, -1 * h/2])
        bottom_right = np.array([w/2, -1 * h/2])

        return np.vstack((bottom_right, top_right, top_left, bottom_left))

    def visualize_minkowski(self, S, ax):
        """Visualize the Minkowski sum S."""
        hull = ConvexHull(S)
        ax.plot(np.append(hull.points[hull.vertices, 0], hull.points[hull.vertices[0], 0]), 
                np.append(hull.points[hull.vertices, 1], hull.points[hull.vertices[0], 1]), color='green')

    #algorithm to compute minkowski sum
    def minkowski_algorithm(self, polygon):
        P = polygon
        Q = -1 * self.rectangle
        
        S = []
        
        for vector1 in P:
            for vector2 in Q:
                S.append(vector1 + vector2)
        
        return np.array(S)
    
    def generate_minkowski_plot(self):
        for polygon in self.polygons:
            self.ax.fill([vertex[0] for vertex in polygon], [vertex[1] for vertex in polygon], alpha=.25, fc='red', ec='black')
            self.visualize_minkowski(self.minkowski_algorithm(polygon), self.ax)

minkowskiPlot = MinkowskiPlot(f, ax, rotation_angle = 0)
minkowskiPlot.generate_minkowski_plot()
plt.show()