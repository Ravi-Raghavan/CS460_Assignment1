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
ax.set_xlim(-1.5, 3.5)
ax.set_ylim(-1.5, 3.5)

class MinkowskiPlot:
    def __init__(self, f, ax):
        #Store figure and axes as instance variables
        self.f = f
        self.ax = ax
        
        #Load Polygon Data
        self.polygons = np.load("rr1133-mmh255/2d_rigid_body.npy", allow_pickle= True)

        #Randomly Generate Starting Configuration for Rectangle
        default_rectangle = self.generate_rectangle() #Generate a rectangle with width w and height h to be centered around the origin
        self.rectangle = self.random_rectangle_configuration(default_rectangle)
        while (self.collides_with_other_polygons(self.rectangle) or self.is_rectangle_on_boundary(self.rectangle)):
            self.rectangle = self.random_rectangle_configuration(default_rectangle)
        
        #Original Rectangle Patch
        self.rectangle_patch = matplotlib.patches.Polygon(self.rectangle, closed=True, facecolor = 'none', edgecolor='r')
        self.ax.add_patch(self.rectangle_patch)

        # Plot centroid of rectangle
        self.body_centroid = self.ax.plot(self.center_point[0],self.center_point[1], marker='o', markersize=3, color="green")
             
    #Rotate and Translate the Default Rectangle by a given angle and translation amount
    def random_rectangle_configuration(self, rectangle):
        self.rotation_angle = np.random.choice(a = np.array([0, 45, 90]))
        self.center_point = np.array([np.random.uniform(0, 2), np.random.uniform(0, 2)])
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

    #Check to see if the rigid body collides with other polygons in our list
    def collides_with_other_polygons(self, rectangle):
        for polygon in self.polygons:
            if (collides(polygon, rectangle)):
                return True
        return False
    
    #Returns true if rectangle is on boundary
    def is_rectangle_on_boundary(self, rectangle):
        for vertex in rectangle:
            if vertex[0] <= 0 or vertex[0] >= 2 or vertex[1] <= 0 or vertex[1] >= 2:
                return True
        
        return False
    
    def visualize_minkowski(self, A, B, S, ax):
        """Visualize polygons A, B, and their Minkowski sum S."""
        # plt.figure(figsize=(8, 8))

        ax.plot(np.append(A[:, 0], A[0, 0]), np.append(A[:, 1], A[0, 1]), label="A", color='blue')
        ax.plot(np.append(B[:, 0], B[0, 0]), np.append(B[:, 1], B[0, 1]), label="B", color='red')

        hull = ConvexHull(S)
        ax.plot(hull.points[hull.vertices, 0], hull.points[hull.vertices, 1], label="Minkowski Sum", color='green')
        ax.plot(np.append(hull.points[hull.vertices, 0], hull.points[hull.vertices[0], 0]), 
                np.append(hull.points[hull.vertices, 1], hull.points[hull.vertices[0], 1]), color='green')
        ax.plot(0, 0, 'bo', label='Origin')

        ax.legend()
        # ax.grid(True)
        # ax.set_aspect('equal')  # Set the aspect ratio of the plot
        # plt.show()

    #compute minkowski sum for set A and B
    def minkowski_sums(self, A, B):
        S = []
        
        for A_vector in A:
            for B_vector in B:
                S.append(A_vector + B_vector)
        
        return np.array(S)

    #compute minkowski difference for set A and B
    def minkowski_difference(self, A, B):
        S = []
        
        for A_vector in A:
            for B_vector in B:
                S.append(A_vector - B_vector)
        
        return np.array(S)
    
    def generate_minkowski_sum_plot(self):
        for polygon in self.polygons:
            self.visualize_minkowski(self.rectangle, polygon, self.minkowski_sums(self.rectangle, polygon), self.ax)

minkowskiPlot = MinkowskiPlot(f, ax)
minkowskiPlot.generate_minkowski_sum_plot()
plt.show()