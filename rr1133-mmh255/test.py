##### test file just to test out certain stuff
import numpy as np


polygon = np.array([[1, 2], [1, 5], [-1, 5], [-1, 2]])
bbox = np.array([np.min(polygon, axis=0), np.max(polygon, axis=0)])

# print(bbox)

F = np.array([1, 2, 3, 4, 5, 1])
print(np.argwhere(F == np.min(F)).flatten())

# T = np.array([[1, 2], [-1, 3], [10,11]])
# print(np.min(T, axis = 0))

# print(np.argmin(F))