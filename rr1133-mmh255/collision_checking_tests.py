import numpy as np 
import matplotlib.pyplot as plt

polygons = np.load("rr1133-mmh255/arm_polygons.npy", allow_pickle= True)
polygon1 = np.vstack((polygons[0], polygons[0][0]))
polygon2 = np.vstack((polygons[1], polygons[1][0]))
polygon3 = np.vstack((polygons[2], polygons[2][0]))

plt.scatter([a[0] for a in polygon1], [a[1] for a in polygon1], c = "blue")
plt.plot([a[0] for a in polygon1], [a[1] for a in polygon1])

plt.scatter([a[0] for a in polygon2], [a[1] for a in polygon2], c = "red")
plt.plot([a[0] for a in polygon2], [a[1] for a in polygon2])

plt.scatter([a[0] for a in polygon3], [a[1] for a in polygon3], c = "black")
plt.plot([a[0] for a in polygon3], [a[1] for a in polygon3])
plt.show()
