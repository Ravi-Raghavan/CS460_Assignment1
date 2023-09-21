import numpy as np 
import matplotlib.pyplot as plt

polygons = np.load("rr1133-mmh255/2d_rigid_body.npy", allow_pickle= True)
fig, ax = plt.subplots(dpi = 100)
ax.set_aspect("equal")

def generate_rectangle():
    center_point = np.array([np.random.uniform(0, 2), np.random.uniform(0, 2)])
    w = 0.2
    h = 0.1
    
    top_left = np.array([-1 * w/2, h/2])
    top_right = np.array([w/2, h/2])
    bottom_left = np.array([-1 * w/2, -1 * h/2])
    bottom_right = np.array([w/2, -1 * h/2])

    angle = np.deg2rad(np.random.uniform(0, 360))
    rotation_matrix = np.array([[np.cos(angle), -1 * np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    
    top_left = (rotation_matrix @ top_left.T).T + center_point
    top_right = (rotation_matrix @ top_right.T).T + center_point

    bottom_left = (rotation_matrix @ bottom_left.T).T + center_point
    bottom_right = (rotation_matrix @ bottom_right.T).T + center_point
    
    return center_point, bottom_right, top_right, top_left, bottom_left

center_point, bottom_right, top_right, top_left, bottom_left = generate_rectangle();
print(bottom_right, top_right, top_left, bottom_left)

print(np.linalg.norm(bottom_right - center_point) * 2)
print(np.linalg.norm(top_right - center_point) * 2)
print(np.linalg.norm(top_left - center_point) * 2)
print(np.linalg.norm(bottom_left - center_point) * 2)

print("------------------------------------------------")

print(np.linalg.norm(bottom_right - top_right))
print(np.linalg.norm(top_right - top_left))
print(np.linalg.norm(top_left - bottom_left))
print(np.linalg.norm(bottom_left - bottom_right))

print("------------------------------------------------")

slope1 = (top_right[1] - bottom_right[1]) / (top_right[0] - bottom_right[0])
slope2 = (bottom_left[1] - bottom_right[1]) / (bottom_left[0] - bottom_right[0])
print(f"{slope1}, {slope2}, {slope1 * slope2}")

slope1 = (top_left[1] - bottom_left[1]) / (top_left[0] - bottom_left[0])
slope2 = (bottom_right[1] - bottom_left[1]) / (bottom_right[0] - bottom_left[0])
print(f"{slope1}, {slope2}, {slope1 * slope2}")

slope1 = (top_right[1] - top_left[1]) / (top_right[0] - top_left[0])
slope2 = (bottom_left[1] - top_left[1]) / (bottom_left[0] - top_left[0])
print(f"{slope1}, {slope2}, {slope1 * slope2}")

slope1 = (top_left[1] - top_right[1]) / (top_left[0] - top_right[0])
slope2 = (bottom_right[1] - top_right[1]) / (bottom_right[0] - top_right[0])
print(f"{slope1}, {slope2}, {slope1 * slope2}")

print("------------------------------------------------")

ax.fill([vertex[0] for vertex in [bottom_right, top_right, top_left, bottom_left]], [vertex[1] for vertex in [bottom_right, top_right, top_left, bottom_left]], alpha=.25, fc='red', ec='blue')

for index in range(len(polygons)):
    polygon = polygons[index]
    ax.fill([vertex[0] for vertex in polygon], [vertex[1] for vertex in polygon], alpha=.25, fc='white', ec='black')
    
plt.show()
