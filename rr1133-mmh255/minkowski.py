import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from collision_checking import *
import platform

#Print python version for sanity check
print(platform.python_version())

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
        self.polygons = np.load("rr1133-mmh255/polygons_scene9.npy", allow_pickle= True)
        
        #Pre-process polygon list to numpy array, just-in-case
        for index, polygon in enumerate(self.polygons):
            polygon = [list(vertex) for vertex in polygon]
            self.polygons[index] = np.array(polygon)

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

    #Visualize Minkowski
    def visualize_minkowski(self, S, ax):
        """Visualize the Minkowski sum S."""
        hull = ConvexHull(S)
        ax.plot(np.append(hull.points[hull.vertices, 0], hull.points[hull.vertices[0], 0]), 
                np.append(hull.points[hull.vertices, 1], hull.points[hull.vertices[0], 1]), color='green')

    #Algorithm to compute minkowski sum
    def minkowski_algorithm(self, polygon):
        P = polygon
        Q = -1 * self.rectangle
        S = []
        
        for vector1 in P:
            for vector2 in Q:
                S.append(vector1 + vector2)
        
        return np.array(S)
    
    #Generating Minkowski Plot
    def generate_minkowski_plot(self):
        for polygon in self.polygons:
            self.ax.fill([vertex[0] for vertex in polygon], [vertex[1] for vertex in polygon], alpha=.25, fc='red', ec='black')
            self.visualize_minkowski(self.minkowski_algorithm(polygon), self.ax)

    
    #Optimized Version of Visualizing Minkowski
    def optimized_visualize_minkowski(self, S, ax):
        ax.plot(np.append(S[:, 0], S[0, 0]), 
                np.append(S[:, 1], S[0, 1]), color='blue')

    
    #Optimized Version of Minkowski Algorithm
    def optimized_minkowski_algorithm(self, polygon):
        #Goal is to compute Minowski Sum of P and Q
        P, Q = polygon, -1 * self.rectangle 
        
        #Get starting pointers for P and Q
        P_pointers, Q_pointers = np.where(P[:, 1] == np.min(P[:, 1]))[0], np.where(Q[:, 1] == np.min(Q[:, 1]))[0]
        P_pointer, Q_pointer = P_pointers[np.argmin(P[P_pointers, 0])], Q_pointers[np.argmin(Q[Q_pointers, 0])]
        P_count, Q_count = 0,0
        
        #Initialize S
        S = []
        previous_P_phi, previous_Q_phi = None, None
                 
        #Iterate
        while P_count < P.shape[0] or Q_count < Q.shape[0]:
            #Sum up the two vertices and store it in S
            Pi, Qj = P[P_pointer % P.shape[0]], Q[Q_pointer % Q.shape[0]]   
            S.append(Pi + Qj)
            
            #Get edges Pi, P_i+1 and edge Qj, Q_j+1
            P_i1, Q_j1 = P[(P_pointer + 1) % P.shape[0]], Q[(Q_pointer + 1) % Q.shape[0]]
            P_edge, Q_edge = P_i1 - Pi, Q_j1 - Qj
            
            #Calculate polar angles. I don't want negative angles so just add 360 degrees. 
            P_phi, Q_phi = np.rad2deg(np.arctan2(P_edge[1], P_edge[0])), np.rad2deg(np.arctan2(Q_edge[1], Q_edge[0]))
            P_phi = P_phi + 360 if P_phi < 0 else P_phi
            Q_phi = Q_phi + 360 if Q_phi < 0 else Q_phi
            
            #Adjust if we have already made a full circle
            P_phi = P_phi + 360 if previous_P_phi != None and P_phi < previous_P_phi else P_phi
            Q_phi = Q_phi + 360 if previous_Q_phi != None and Q_phi < previous_Q_phi else Q_phi
                        
            #Adjust Pointer
            if P_phi < Q_phi:
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
            
            previous_P_phi, previous_Q_phi = P_phi, Q_phi
                         
        return np.array(S)
        
    #Optimized verison of generation of Minkowski Plot
    def optimized_generate_minkowski_plot(self):
        for polygon in self.polygons:
            self.ax.fill([vertex[0] for vertex in polygon], [vertex[1] for vertex in polygon], alpha=.25, fc='red', ec='black')
            self.optimized_visualize_minkowski(self.optimized_minkowski_algorithm(polygon), self.ax)

minkowskiPlot = MinkowskiPlot(f, ax, rotation_angle = 0)
minkowskiPlot.optimized_generate_minkowski_plot()
minkowskiPlot.generate_minkowski_plot()
plt.show()