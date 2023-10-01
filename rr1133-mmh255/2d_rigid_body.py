import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from collision_checking import *
from scipy.spatial import ConvexHull
import argparse

#Print Matplotlib version for sanity
print(f"matplotlib version: {matplotlib.__version__}")

# Python Class For 2D Rigid Body
class RigidBody:
    def __init__(self, f, ax, file):
        #Store figure and axes as instance variables
        self.f = f
        self.ax = ax
        
        #set axis limits for x and y axes
        ax.set_xlim(0, 2)
        ax.set_ylim(0, 2)
        
        #set title for axis
        ax.set_title('2D Rigid Body Simulation', fontsize = 12)
        
        #Load Polygon Data
        self.polygons = np.load(file, allow_pickle= True)

        #Plot polygons from the npy file
        for index in range(len(self.polygons)):
            self.ax.fill([vertex[0] for vertex in self.polygons[index]], [vertex[1] for vertex in self.polygons[index]], alpha=.25, fc='white', ec='black')

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

        # Connect keyboard events to event handlers
        self.f.canvas.mpl_connect('key_press_event', self.keyboard_event_handler)

    #Returns true if rectangle is on boundary
    def is_rectangle_on_boundary(self, rectangle):
        for vertex in rectangle:
            if vertex[0] <= 0 or vertex[0] >= 2 or vertex[1] <= 0 or vertex[1] >= 2:
                return True
        return False
    
    #Rotate rectangle about center
    def rotate_about_center(self, event):
        #Change angle of rotation
        delta_rotation_angle = 10 if event == "up" else -10
        new_rotation_angle = self.rotation_angle + delta_rotation_angle
        new_rotation_angle = new_rotation_angle + 360 if new_rotation_angle < 0 else new_rotation_angle
        new_rotation_angle = new_rotation_angle % 360
        
        #Rotate the Rectangle
        default_rectangle = self.generate_rectangle()
        angle = np.deg2rad(new_rotation_angle)
        rotation_matrix = np.array([[np.cos(angle), -1 * np.sin(angle)], [np.sin(angle), np.cos(angle)]])  
        new_rectangle = (rotation_matrix @ default_rectangle.T).T + self.center_point
        
        if not (self.collides_with_other_polygons(new_rectangle) or self.is_rectangle_on_boundary(new_rectangle)):
            self.rectangle = new_rectangle
            self.rotation_angle = new_rotation_angle
    
    #Rewriting the rotation algorithm
    def rotate_about_center_2(self, event):
        #Change angle of rotation
        delta_rotation_angle = 10 if event == "up" else -10
        
        #Set up new rotation angle
        new_rotation_angle = self.rotation_angle + delta_rotation_angle
        new_rotation_angle = new_rotation_angle + 360 if new_rotation_angle < 0 else new_rotation_angle
        new_rotation_angle = new_rotation_angle % 360
                
        #Modify the Rectangle
        modified_rectangle = self.rectangle - self.center_point
        modified_rectangle = np.hstack((modified_rectangle, np.ones(shape = (modified_rectangle.shape[0], 1)))).T

        #Rotate the Rectangle
        angle = np.deg2rad(delta_rotation_angle)
        rotation_matrix = np.array([[np.cos(angle), -1 * np.sin(angle), self.center_point[0]], [np.sin(angle), np.cos(angle), self.center_point[1]], [0, 0, 1]])
        new_rectangle = (rotation_matrix @ modified_rectangle).T 
        new_rectangle = new_rectangle[:, :-1]
        
        if not (self.collides_with_other_polygons(new_rectangle) or self.is_rectangle_on_boundary(new_rectangle)):
            self.rectangle = new_rectangle
            self.rotation_angle = new_rotation_angle

    def translate_rectangle(self, event):
        #r: total translation across x and y
        r = 0.05
        angle_of_rotation = self.rotation_angle if event == "right" else self.rotation_angle - 180
        delta = np.array([r * np.cos(np.deg2rad(angle_of_rotation)), r * np.sin(np.deg2rad(angle_of_rotation)), 1])
        
        #Form the Translation Matrix
        translation_matrix = np.identity(3)
        translation_matrix[:, -1] = delta
        
        #Modify the Center of the Rectangle using the translation matrix
        modified_center_point = np.vstack((self.center_point.reshape(-1, 1), np.array([[1]])))
        modified_center_point = translation_matrix @ modified_center_point
        modified_center_point = modified_center_point.flatten()[:-1]
        
        #Modify the Rectangle Vertices using the translation matrix
        modified_rectangle = np.hstack((self.rectangle, np.ones(shape = (self.rectangle.shape[0], 1)))).T
        modified_rectangle = (translation_matrix @ modified_rectangle).T
        modified_rectangle = modified_rectangle[:, :-1]
        
        if not (self.collides_with_other_polygons(modified_rectangle) or self.is_rectangle_on_boundary(modified_rectangle)):
            self.center_point = modified_center_point
            self.rectangle = modified_rectangle
        
        self.body_centroid[0].set_data([self.center_point[0], self.center_point[1]])
             
    #Rotate and Translate the Default Rectangle by a given angle and translation amount
    def random_rectangle_configuration(self, rectangle):
        self.rotation_angle = np.random.uniform(0, 360)
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
            if (collides_optimized(polygon, np.vstack((rectangle, rectangle[0])))):
                return True
        return False
    
    # Event handler to change the rotation angle
    def keyboard_event_handler(self, event):
        if event.key == "up" or event.key == 'down':
            self.rotate_about_center_2(event.key)
            self.rectangle_patch.set_xy(self.rectangle)
        elif event.key == 'right' or event.key == 'left':
            self.translate_rectangle(event.key)
            self.rectangle_patch.set_xy(self.rectangle)
        
        self.f.canvas.draw()


class MinkowskiPlot:
    def __init__(self, f, ax, rotation_angle, file):
        #Store figure and axes as instance variables
        self.f = f
        self.ax = ax
        self.rotation_angle = rotation_angle
        
        #set axis limits for x and y axes
        ax.set_xlim(-0.5, 2.5)
        ax.set_ylim(-0.5, 2.5)
        
        #Load Polygon Data
        self.polygons = np.load(file, allow_pickle= True)
        
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
        P = polygon[:-1]
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
        P, Q = polygon[:-1], -1 * self.rectangle 
        
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

# Create an argument parser
parser = argparse.ArgumentParser(description='2D Rigid Body Python Program')

# Add a string argument for "function"
parser.add_argument('--function', type=str, help='function parameter')

# Parse the command-line arguments
args = parser.parse_args()

# Access the "function" argument
function_parameter = args.function

if function_parameter == "Rigid":
    #Initialize figures and axes
    f,ax = plt.subplots(dpi = 100)
    ax.set_aspect("equal")
    
    file = "rr1133-mmh255/2d_rigid_body.npy"
    rigidBody = RigidBody(f, ax, file)
    plt.show()
elif function_parameter == "Minkowski":
    #Initialize figures and axes
    f,ax = plt.subplots(dpi = 100)
    for scene_number in range(1, 9):
        file = f"rr1133-mmh255/scene.npy" if scene_number == 1 else f"rr1133-mmh255/scene{scene_number}.npy"        
        for rotation_angle in [0, 45, 90]:
            ax.set_aspect("equal")
            minkowskiPlot = MinkowskiPlot(f, ax, rotation_angle = rotation_angle, file = file)
            minkowskiPlot.optimized_generate_minkowski_plot()
            minkowskiPlot.generate_minkowski_plot()
            plt.savefig(f"Problem3_minkowski{scene_number}_{rotation_angle}.jpg")
            ax.clear()

    plt.close(f)