import numpy as np 
import matplotlib.pyplot as plt

polygons = np.load("rr1133-mmh255/arm_polygons.npy", allow_pickle= True)

for polygon in polygons:
    print("POLYGON")
    print(polygon)
    print(polygon.shape[0])
    print("------------")
    

#given a polygon as numpy array, get its edges as an array
def get_edges(polygon):
    V = polygon.shape[0]
    edges = [[polygon[i], polygon[(i + 1) % V]] for i in range(V)]
    return edges

def determinant(A):
    return (A[1, 1] * A[0, 0]) - (A[0, 1] * A[1, 0])

#Given two edges, determine if they intersect using parameterization and Cramer's Rule
def edge_intersect(edge1, edge2):
    ## Linear Segments can be written using s and t parameters. Edge 1 will be expressed using s parameter and Edge 2 will be expressed using t parameter
    
    #Edge 1 x and y values
    x1 = edge1[0][0]
    x2 = edge1[1][0]
    y1 = edge1[0][1]
    y2 = edge1[1][1]
    
    #Edge 2 x and y values
    x3 = edge2[0][0]
    x4 = edge2[1][0]
    y3 = edge2[0][1]
    y4 = edge2[1][1]
    
    #Determinants
    Dx = determinant(np.array([[x3 - x1, x3 - x4], [y3 - y1, y3 - y4]]))
    Dy = determinant(np.array([[x2 - x1, x3 - x1], [y2 - y1, y3 - y1]]))
    D = determinant(np.array([[x2 - x1, x3 - x4], [y2 - y1, y3 - y4]]))
    
    if D == 0:
        return Dx == 0 and Dy == 0    
    
    s = Dx / D
    t = Dy / D
    
    return 0 <= s and s <= 1 and 0 <= t and t <= 1    

## EXHAUSTIVE APPROACH FOR DETECTING POLYGON COLLISION##
def collides(poly1, poly2):
    poly1_edges = get_edges(poly1)
    poly2_edges = get_edges(poly2)
    
    for edge1 in poly1_edges:
        for edge2 in poly2_edges:
            if (edge_intersect(edge1, edge2)):
                return True

    return False

## Optimized Approach for Detecting Polygon Collision
def collides_optimized(poly1, poly2):
    pass


def plot(polygons):
    pass