import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from collision_checking import *

print(f"matplotlib version: {matplotlib.__version__}")
polygons = np.load("rr1133-mmh255/2d_rigid_body.npy", allow_pickle= True)

rotation_angle = None
center_point = None

f,ax = plt.subplots(dpi = 100)
ax.set_aspect("equal")
ax.set_xlim(-0.01, 2.01)
ax.set_ylim(-0.01, 2.01)


for index in range(len(polygons)):
    polygon = polygons[index]
    ax.fill([vertex[0] for vertex in polygon], [vertex[1] for vertex in polygon], alpha=.25, fc='white', ec='black')

def is_rectangle_on_boundary(rectangle):
    for vertex in rectangle:
        if (vertex[0] <= 0 or vertex[0] >= 2):
            return True
        
        if vertex[1] <= 0 or vertex[1] >= 2:
            return True
    
    return False
    
def rotate_about_center(rectangle, event):
    global rotation_angle, center_point
    original_rectangle = rectangle
    original_rotation_angle = rotation_angle
    original_center_point = center_point 
    
    print("------------------ROTATION-----------------------")
    print(f"Start, Rotation Angle: {rotation_angle}, Center Point: {center_point}, bottom_right: {rectangle[0]}, bottom_left: {rectangle[-1]}")
    
    delta_rotation_angle = 10 if event == "up" else -10
    rotation_angle += delta_rotation_angle
    rotation_angle = rotation_angle % 360
    
    rectangle = generate_rectangle()
    angle = np.deg2rad(rotation_angle)
    rotation_matrix = np.array([[np.cos(angle), -1 * np.sin(angle)], [np.sin(angle), np.cos(angle)]])  
    rectangle = (rotation_matrix @ rectangle.T).T + center_point
    
    if (collides_with_other_polygons(rectangle) or is_rectangle_on_boundary(rectangle)):
        rectangle = original_rectangle
        rotation_angle = original_rotation_angle
        center_point = original_center_point
    
    print(f"End, Rotation Angle: {rotation_angle}, Center Point: {center_point}, bottom_right: {rectangle[0]}, bottom_left: {rectangle[-1]}")
    print("-------------------------------------------------")

    return rectangle

def translate_rectangle(rectangle, event):
    global rotation_angle, center_point
    original_rectangle = rectangle
    original_rotation_angle = rotation_angle
    original_center_point = center_point 
    
    print("------------------TRANSLATION-----------------------")
    print(f"Start, Rotation Angle: {rotation_angle}, Center Point: {center_point}, bottom_right: {rectangle[0]}, bottom_left: {rectangle[-1]}")
    
    r = 0.05
    delta_center_point = None 
    
    if event == "right":
        delta_cx = r * np.cos(np.deg2rad(rotation_angle))
        delta_cy = r * np.sin(np.deg2rad(rotation_angle))
        delta_center_point = np.array([delta_cx, delta_cy])
    
    elif event == "left":
        delta_cx = r * np.cos(np.deg2rad(rotation_angle - 180))
        delta_cy = r * np.sin(np.deg2rad(rotation_angle - 180))
        delta_center_point = np.array([delta_cx, delta_cy])

    center_point += delta_center_point
    rectangle = rectangle + delta_center_point
    
    if (collides_with_other_polygons(rectangle) or is_rectangle_on_boundary(rectangle)):
        rectangle = original_rectangle
        rotation_angle = original_rotation_angle
        center_point = original_center_point
    
    print(f"End, Rotation Angle: {rotation_angle}, Center Point: {center_point}, bottom_right: {rectangle[0]}, bottom_left: {rectangle[-1]}")
    print("-------------------------------------------------")
    
    return rectangle
     
#Rotate and Translate the Rectangle by a given angle and translation amount
def randomly_rotate_and_translate_rectangle(rectangle):
    global rotation_angle, center_point
    rotation_angle = np.random.uniform(0, 360)
    
    center_point = np.array([np.random.uniform(0, 2), np.random.uniform(0, 2)])
    angle = np.deg2rad(rotation_angle)
    
    rotation_matrix = np.array([[np.cos(angle), -1 * np.sin(angle)], [np.sin(angle), np.cos(angle)]])    
    return (rotation_matrix @ rectangle.T).T + center_point

#Generate Rectangle that is centered around the origin
def generate_rectangle():
    w = 0.2
    h = 0.1
    
    top_left = np.array([-1 * w/2, h/2])
    top_right = np.array([w/2, h/2])
    bottom_left = np.array([-1 * w/2, -1 * h/2])
    bottom_right = np.array([w/2, -1 * h/2])

    return np.vstack((bottom_right, top_right, top_left, bottom_left))

def collides_with_other_polygons(rectangle):
    for polygon in polygons:
        if (collides(polygon, rectangle)):
            return True
    return False

default_rectangle = generate_rectangle()
rectangle = randomly_rotate_and_translate_rectangle(default_rectangle)
while (collides_with_other_polygons(rectangle)):
    rectangle = randomly_rotate_and_translate_rectangle(default_rectangle)

coords = center_point #center of rectangle 

#Original Rectangle Patch
rec0 = matplotlib.patches.Polygon(rectangle, closed=True, facecolor = 'none', edgecolor='r')
patch = ax.add_patch(rec0)

# Rectangles center
plt.plot(coords[0],coords[1], marker='o', markersize=3, color="green")

# Event handler to change the rotation angle
def change_rotation_or_translation(event):
    global rectangle, rotation_angle
    if event.key == "up" or event.key == 'down':
        rectangle = rotate_about_center(rectangle, event.key)
        print(f"RECTANGLE: {rectangle}")
        rec0.set_xy(rectangle)
    elif event.key == 'right' or event.key == 'left':
        rectangle = translate_rectangle(rectangle, event.key)
        print(f"RECTANGLE: {rectangle}")
        rec0.set_xy(rectangle)
    
    f.canvas.draw()

# Connect keyboard events to event handlers
f.canvas.mpl_connect('key_press_event', change_rotation_or_translation)
plt.show()