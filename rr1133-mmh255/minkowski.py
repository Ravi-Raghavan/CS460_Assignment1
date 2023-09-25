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
    
    def optimized_visualize_minkowski(self, S, ax):
        ax.plot(np.append(S[:, 0], S[0, 0]), 
                np.append(S[:, 1], S[0, 1]), color='green')
    
    #Optimized Version of Minkowski Algorithm
    def optimized_minkowski_algorithm(self, polygon):
        P = polygon
        Q = -1 * self.rectangle
        
        P_pointer, Q_pointer = np.argmin(P[:, 1]), np.argmin(Q[:, 1])
        P_count, Q_count = 0,0
        
        S = []
        
        while P_count < P.shape[0] or Q_count < Q.shape[0]:
            Pi, Qj = P[P_pointer % P.shape[0]], Q[Q_pointer % Q.shape[0]]            
            S.append(Pi + Qj)
            
            P_i1, Q_j1 = P[(P_pointer + 1) % P.shape[0]], Q[(Q_pointer + 1) % Q.shape[0]]
            P_edge, Q_edge = P_i1 - Pi, Q_j1 - Qj
            P_phi, Q_phi = np.rad2deg(np.arctan2(P_edge[1], P_edge[0])), np.rad2deg(np.arctan2(Q_edge[1], Q_edge[0]))
            
            P_phi = P_phi + 360 if P_phi < 0 else P_phi
            Q_phi = Q_phi + 360 if Q_phi < 0 else Q_phi

            if (P_phi < Q_phi):
                P_pointer += 1
                P_count += 1
            elif P_phi > Q_phi:
                Q_pointer += 1
                Q_count += 1
            else:
                P_pointer += 1
                P_count += 1
                Q_pointer += 1
                Q_count += 1
                
            print(Pi, Qj, P_phi, Q_phi)
            
            if (P_count >= 10 or Q_count >= 10):
                break
         
                    
        return np.array(S)
    
    def generate_minkowski_plot(self):
        for polygon in self.polygons:
            self.ax.fill([vertex[0] for vertex in polygon], [vertex[1] for vertex in polygon], alpha=.25, fc='red', ec='black')
            self.visualize_minkowski(self.minkowski_algorithm(polygon), self.ax)
    
    def optimized_generate_minkowski_plot(self):
        for polygon in self.polygons[0:1]:
            self.ax.fill([vertex[0] for vertex in polygon], [vertex[1] for vertex in polygon], alpha=.25, fc='red', ec='black')
            self.optimized_visualize_minkowski(self.optimized_minkowski_algorithm(polygon), self.ax)

minkowskiPlot = MinkowskiPlot(f, ax, rotation_angle = 0)
minkowskiPlot.optimized_generate_minkowski_plot()
plt.show()