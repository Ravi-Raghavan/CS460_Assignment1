import numpy as np 
import matplotlib.pyplot as plt

polygons = np.load("rr1133-mmh255/collision_checking_polygons.npy", allow_pickle= True)
polygon1 = polygons[0]
polygon2 = polygons[1]
polygon3 = polygons[2]

plt.scatter([a[0] for a in polygon1], [a[1] for a in polygon1], c = "blue")
plt.plot([a[0] for a in polygon1], [a[1] for a in polygon1])

plt.scatter([a[0] for a in polygon2], [a[1] for a in polygon2], c = "red")
plt.plot([a[0] for a in polygon2], [a[1] for a in polygon2])

plt.scatter([a[0] for a in polygon3], [a[1] for a in polygon3], c = "black")
plt.plot([a[0] for a in polygon3], [a[1] for a in polygon3])
plt.show()
