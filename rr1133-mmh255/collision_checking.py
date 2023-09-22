import numpy as np 
import matplotlib.pyplot as plt

#given a polygon as numpy array, get its edges as an array
def get_edges(polygon):
    V = polygon.shape[0]
    edges = [[polygon[i], polygon[i + 1]] for i in range(V - 1)]
    return edges

#calculate determinant of 2 x 2 matrix
def determinant(A):
    return (A[1, 1] * A[0, 0]) - (A[0, 1] * A[1, 0])

#calculate Euclidean Distance between two points
def euclidean_distance(point1, point2):
    return np.sqrt((point2[1] - point1[1]) ** 2 + (point2[0] - point1[0]) ** 2)

#determine if edge is on point
def point_on_edge(edge, point, epsilon = 1e-50):
    return (euclidean_distance(edge[0], point) + euclidean_distance(edge[1], point) - euclidean_distance(edge[0], edge[1])) < epsilon

#Given two polygons, check to see if a point is contained within a polygon
def point_contained(poly, point):
    center = poly[0]
    angles = []
    special_case_phi = False #Case where our angles go from 4th quadrant to 1st quadrant
    
    for vertex in poly[1: -1]:
        relative_point = vertex - center
        phi = np.rad2deg(np.arctan2(relative_point[1], relative_point[0]))
        
        #keep value of phi between 0 and 360
        phi = phi + 360 if phi < 0 else phi
        
        #If values start to decrease, we know we are going from 4th quadrant to 1st quadrant. In this case, simply add 360
        if (len(angles) >= 1 and phi < angles[-1]):
            phi += 360  
            special_case_phi = True
            
        angles.append(phi)
    
    relative_point = point - center
    point_phi = np.rad2deg(np.arctan2(relative_point[1], relative_point[0]))
    point_phi = point_phi + 360 if point_phi < 0 or special_case_phi else point_phi
    
    index = np.searchsorted(angles, point_phi)
    if (index == len(angles) or (index == 0 and point_phi < angles[0])):
        return False
    
    #If angle matches exactly, shortcut check
    if (point_phi == angles[index]):
        return euclidean_distance(center, point) <= euclidean_distance(center, poly[index + 1])
    
    ray = [center, point]
    edge = [poly[index], poly[index + 1]]
            
    return (not edge_intersect(ray, edge)) or point_on_edge(edge, point)
        

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