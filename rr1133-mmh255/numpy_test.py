import numpy as np

A = np.array([[1, 2], [3, 4]])
B = np.array([10, 20])
print(A + B)

print(np.rad2deg(np.arctan2(1.732, 1)))
print(np.searchsorted([1,2,3,4,5], 5))


print(B.reshape(1,-1))

C = np.vstack((B.reshape(1,-1).T, np.array([[1]])))
print(C.flatten())


B = np.array([1, 2, 3, 4, 5])
print(B.reshape((-1, 1)))


M = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
D = np.array([100, 200, 300])
T = np.array([[1], [2], [3]])
M[:, -1] = D
print(M)

print(T.flatten()[: -1])

C = np.array([1, 1, 2, 3, 4])
print(np.argmin(C))
print(np.where(C == 1))


D = np.array([[1, 2], [3, -1], [5, -1]])
print(D[:, 1].flatten() == np.min(D[:, 1]))
print(np.where(D[:, 1].flatten() == np.min(D[:, 1].flatten())))