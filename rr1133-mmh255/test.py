import matplotlib.pyplot as plt
import matplotlib

print(f"matplotlib version: {matplotlib.__version__}")

f,ax = plt.subplots(dpi = 100)
ax.set_aspect("equal")
coords = [2.5,2] #center of rectangle 

#Original Rectangle Patch
rec0 = matplotlib.patches.Rectangle((1,1),3,2,linewidth=1,edgecolor='r',facecolor='none')
ax.add_patch(rec0)

#Rotated rectangle patch
rect1 = matplotlib.patches.Rectangle((1,1),3,2,linewidth=1,edgecolor='b',facecolor='none',angle = 20, rotation_point = 'center')
ax.add_patch(rect1)

# Rectangles center
plt.plot(coords[0],coords[1], marker='o', markersize=3, color="green")

plt.grid(True)
ax.set_xlim(0,6)
ax.set_ylim(-1,4)

plt.show()
