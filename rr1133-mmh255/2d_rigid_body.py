import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from collision_checking import *

#Print Matplotlib version for sanity
print(f"matplotlib version: {matplotlib.__version__}")

#Initialize figures and axes
f,ax = plt.subplots(dpi = 100)
ax.set_aspect("equal")

#set axis limits for x and y axes
ax.set_xlim(-0.01, 2.01)
ax.set_ylim(-0.01, 2.01)

# Python Class For 2D Rigid Body
class RigidBody:
    def __init__(self, f, ax):
        #Store figure and axes as instance variables
        self.f = f
        self.ax = ax
        
        #Load Polygon Data
        self.polygons = np.load("rr1133-mmh255/2d_rigid_body.npy", allow_pickle= True)

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
        print("------------------ROTATION-----------------------")
        print(f"Start, Rotation Angle: {self.rotation_angle}, Center Point: {self.center_point}, bottom_right: {self.rectangle[0]}, bottom_left: {self.rectangle[-1]}")
        
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
        
        print(f"End, Rotation Angle: {self.rotation_angle}, Center Point: {self.center_point}, bottom_right: {self.rectangle[0]}, bottom_left: {self.rectangle[-1]}")
        print("-------------------------------------------------")

    def translate_rectangle(self, event):
        print("------------------TRANSLATION-----------------------")
        print(f"Start, Rotation Angle: {self.rotation_angle}, Center Point: {self.center_point}, bottom_right: {self.rectangle[0]}, bottom_left: {self.rectangle[-1]}")
        
        #r: total translation across x and y
        r = 0.05
        delta_center_point = None 
        
        #right means g
        if event == "right":
            delta_cx = r * np.cos(np.deg2rad(self.rotation_angle))
            delta_cy = r * np.sin(np.deg2rad(self.rotation_angle))
            delta_center_point = np.array([delta_cx, delta_cy])
        elif event == "left":
            delta_cx = r * np.cos(np.deg2rad(self.rotation_angle - 180))
            delta_cy = r * np.sin(np.deg2rad(self.rotation_angle - 180))
            delta_center_point = np.array([delta_cx, delta_cy])

        self.center_point += delta_center_point
        self.rectangle += delta_center_point
        
        if (self.collides_with_other_polygons(self.rectangle) or self.is_rectangle_on_boundary(self.rectangle)):
            self.center_point -= delta_center_point
            self.rectangle -= delta_center_point
        
        self.body_centroid[0].set_data([self.center_point[0], self.center_point[1]])
        print(f"End, Rotation Angle: {self.rotation_angle}, Center Point: {self.center_point}, bottom_right: {self.rectangle[0]}, bottom_left: {self.rectangle[-1]}")
        print("-------------------------------------------------")
             
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
            if (collides(polygon, rectangle)):
                return True
        return False
    
    # Event handler to change the rotation angle
    def keyboard_event_handler(self, event):
        if event.key == "up" or event.key == 'down':
            self.rotate_about_center(event.key)
            print(f"RECTANGLE: {self.rectangle}")
            self.rectangle_patch.set_xy(self.rectangle)
        elif event.key == 'right' or event.key == 'left':
            self.translate_rectangle(event.key)
            print(f"RECTANGLE: {self.rectangle}")
            self.rectangle_patch.set_xy(self.rectangle)
        
        self.f.canvas.draw()


rigidBody = RigidBody(f, ax)
plt.show()