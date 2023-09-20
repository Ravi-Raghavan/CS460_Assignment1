import numpy as np 
import matplotlib.pyplot as plt

polygons = np.load("rr1133-mmh255/arm_polygons.npy", allow_pickle= True)

#given a polygon as numpy array, get its edges as an array
def get_edges(polygon):
    V = polygon.shape[0]
    edges = [[polygon[i], polygon[(i + 1) % V]] for i in range(V)]
    return edges

#calculate determinant of 2 x 2 matrix
def determinant(A):
    return (A[1, 1] * A[0, 0]) - (A[0, 1] * A[1, 0])

#Given two polygons, check to see if a point is contained within a polygon
def point_contained(poly, point):
    center = poly[0]
    angles = []
    
    for vertex in poly[1: ]:
        relative_point = vertex - center
        phi = np.rad2deg(np.arctan2(relative_point[1], relative_point[0]))
        phi = phi + 360 if phi < 0 else phi            
        angles.append(phi)
    
    relative_point = point - center
    point_phi = np.rad2deg(np.arctan2(relative_point[1], relative_point[0]))
    point_phi = point_phi + 360 if point_phi < 0 else point_phi
    
    index = np.searchsorted(angles, point_phi)
    if ((index > len(angles) and point_phi > angles[-1]) or (index == 0 and point_phi < angles[0])):
        return False
    
    ray = [center, point]
    edge = [poly[0], poly[1]] if index == 0 else [poly[index - 1], poly[index]]    
    return edge_intersect(ray, edge)
        

#Given two edges, determine if they intersect using parameterization and Cramer's Rule
def edge_intersect(edge1, edge2):
    ## Linear Segments can be written using s and t parameters. Edge 1 will be expressed using s parameter and Edge 2 will be expressed using t parameter
    
    #Edge 1 x and y values
    x1, x2, y1, y2 = edge1[0][0], edge1[1][0], edge1[0][1], edge1[1][1]
    
    #Edge 2 x and y values
    x3, x4, y3, y4 = edge2[0][0], edge2[1][0], edge2[0][1], edge2[1][1]
    
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
    #Check to see if each vertex in poly1 is contained in poly2
    for vertex in poly1:
        if (point_contained(poly2, vertex)):
            return True
    
    #Check for Edge Intersection between both polygons
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


#Test Code to test functions I wrote
polygon1 = polygons[0]
polygon2 = polygons[1]
polygon3 = polygons[2]

print(polygon1)
print("-------------------------------")
print(polygon2)
print("-------------------------------")
print(polygon3)
print("-------------------------------")

print(collides(polygon1, polygon2))
print(collides(polygon1, polygon3))
print(collides(polygon3, polygon2))
