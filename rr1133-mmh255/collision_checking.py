import numpy as np 
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

#given a polygon as numpy array, get its edges as an array
#CONFIRMED WORKS
#Expected Format of Input: "polygon" is an (n + 1) x 2 array. There are n vertices in the polygon but the 1st vertex has to be repeated at the end to indicate polygon is closed
#Returns: list of form [edge 1, edge 2, ...., edge n] where edge i is of the form [vertex i, vertex i + 1]
def get_edges(polygon):
    V = polygon.shape[0]
    
    #subtract 1 because first and last vertex is repeated
    edges = [[polygon[i], polygon[i + 1]] for i in range(V - 1)]
    return edges

#calculate determinant of 2 x 2 Matrix A
#CONFIRMED WORKS
def determinant(A):
    return (A[1, 1] * A[0, 0]) - (A[0, 1] * A[1, 0])

#calculate Euclidean Distance between two points
#point1 is of form [x1, y1] and point2 is of form [x2, y2]
#CONFIRMED WORKS
def euclidean_distance(point1, point2):
    return np.sqrt((point2[1] - point1[1]) ** 2 + (point2[0] - point1[0]) ** 2)

#determine if edge is on point
#edge is of form [vertex i, vertex i + 1]
#Each vertex is represented as a point of the form [x1, y1]
#point is of form [x1, y1]
#CONFIRMED WORKS
def point_on_edge(edge, point):
    point1, point2 = edge[0], edge[1]
    m = (point2[1] - point1[1]) / (point2[0] - point1[0])
    
    #x range on the edge
    min_x, max_x = min(point1[0], point2[0]), max(point1[0], point2[0])
    
    #y range on edge    
    min_y, max_y = min(point1[1], point2[1]), max(point1[1], point2[1])
    
    b = point2[1] - (m * point2[0])
    
    inRange = min_x <= point[0] and point[0] <= max_x and min_y <= point[1] and point[1] <= max_y
    onEdge = (point[1] == b + m * point[0])
    return inRange and onEdge

#Given two polygons, check to see if a point is contained within a polygon
#"polygon" is an (n + 1) x 2 array. There are n vertices in the polygon but the 1st vertex has to be repeated at the end to indicate polygon is closed
#"point" is of form [x1, y1]
#CONFIRMED WORKS
def point_contained(poly, point):
    center = poly[0] #center of polar system we are using for reference
    angles = [] #store list of angles we calculate to divide the polygon into sections
    passed_full_circle = False #keep track if we pass full circle
    
    #iterate through all the vertices aside from 'center'
    for vertex in poly[1: -1]:
        relative_point = vertex - center
        phi = np.rad2deg(np.arctan2(relative_point[1], relative_point[0]))
        
        #keep value of phi between 0 and 360
        phi = phi + 360 if phi < 0 else phi
        
        #If values start to decrease, we have passed the 0 degree mark and made a full circle. In this case, simply add 360
        if (len(angles) >= 1 and phi < angles[-1]):
            phi += 360  
            passed_full_circle = True
            
        angles.append(phi)
    
    #Calculate vector from center point to the point of interest. Also calculate the angle this makes with respect to the horizontal axis emitted from our center point
    relative_point = point - center
    point_phi = np.rad2deg(np.arctan2(relative_point[1], relative_point[0]))
    
    #Adjust the angle of point_phi accordingly
    point_phi = point_phi + 360 if point_phi < 0 else point_phi #First get it within the [0, 360] range    
    point_phi = point_phi + 360 if passed_full_circle and (point_phi < angles[0] or point_phi > angles[-1]) else point_phi #Now adjust for the fact that we may have gone past the origin in our circle
    
    #Apply a binary search
    index = np.searchsorted(angles, point_phi)
    if (index == len(angles) or (index == 0 and point_phi < angles[0])):
        return False
    
    #If angle matches exactly, shortcut check
    if (point_phi == angles[index]):
        return euclidean_distance(center, point) <= euclidean_distance(center, poly[index + 1])
    
    #If not, we must do a ray, edge check
    ray = [center, point]
    edge = [poly[index], poly[index + 1]]
            
    return (not edge_intersect(ray, edge)) or point_on_edge(edge, point)
        

#Given two edges, determine if they intersect using parameterization and Cramer's Rule
#CONFIRMED WORKS
def edge_intersect(edge1, edge2):
    ## Linear Segments can be written using s and t parameters. Edge 1 will be expressed using s parameter and Edge 2 will be expressed using t parameter
    
    #Edge 1 x and y values. Edge 1 connects [x1, y1] to [x2, y2]
    x1, x2, y1, y2 = edge1[0][0], edge1[1][0], edge1[0][1], edge1[1][1]
    
    #Edge 2 x and y values. Edge 2 connects [x3 y3] to [x4, y4]
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
#CONFIRMED WORKS
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


#Check for collision between two bounding boxes
#CONFIRMED WORKS
def check_collision(bbox1, bbox2):
    return not (bbox1[1][0] < bbox2[0][0] or 
                bbox1[0][0] > bbox2[1][0] or 
                bbox1[1][1] < bbox2[0][1] or 
                bbox1[0][1] > bbox2[1][1])
    

#returns True if a collision has been detected, else returns False
#Note: both polygons are an (n + 1) x 2 array. There are n vertices in the polygon but the 1st vertex has to be repeated at the end to indicate polygon is closed
#CONFIRMED WORKS
def SAT(poly1, poly2):
    poly1_edges = get_edges(poly1)
    poly2_edges = get_edges(poly2)
    edges = poly1_edges + poly2_edges
    
    for edge in edges:
        edge_vector = edge[1] - edge[0]
        normal_vector = np.array([-1 * edge_vector[1], edge_vector[0]])
        normal_vector /= np.linalg.norm(normal_vector)
        
        poly1_projections = []
        poly2_projections = []
        
        #Compute Projections
        for vertex in poly1:
            projection = np.dot(vertex, normal_vector) * normal_vector
            poly1_projections.append(projection[0])
        
        for vertex in poly2:
            projection = np.dot(vertex, normal_vector) * normal_vector
            poly2_projections.append(projection[0])
        
        poly1_projections = np.sort(np.array(poly1_projections))
        poly2_projections = np.sort(np.array(poly2_projections))
        
        if (poly2_projections[-1] < poly1_projections[0] or poly2_projections[0] > poly1_projections[-1]):
            return False
    
    return True

## Optimized Approach for Detecting Polygon Collision
#Note: both polygons are an (n + 1) x 2 array. There are n vertices in the polygon but the 1st vertex has to be repeated at the end to indicate polygon is closed
#CONFIRMED WORKS
def collides_optimized(poly1, poly2):
    bounding_boxes = [np.array([np.min(polygon, axis=0), np.max(polygon, axis=0)]) for polygon in [poly1, poly2]]
    if (check_collision(bounding_boxes[0], bounding_boxes[1])):
        return SAT(poly1, poly2)

## Given polygons in a scene, use the collision detection algorithm to color them if they collide. If they don't collide with anything, don't color the polygon. 
#CONFIRMED WORKS
def plot(polys, output_file_name = "Problem2_scene.jpg", display_plot = False, save_plot = False, print_diagnostics = False, collision_detection_function = "Minkowski"):
    #Code to be used for plotting figures as assignment requests
    if save_plot or display_plot:
        fig, ax = plt.subplots(dpi = 100)
        ax.set_aspect("equal")

    #for each polygon check to see if collides with other polygons
    for index1 in range(len(polys)):
        polygon = polys[index1]
        collision_polygon = False
        
        # check the other polygons to see if it collides
        for index2 in range(len(polys)):
            
            #don't compare same polygon with itself
            if (index1 == index2):
                continue
            
            #call collision detection algorithm
            if collision_detection_function == "Minkowski":
                if (collides_optimized_alternative(polygon, polys[index2])):
                    collision_polygon = True
                    break
            elif collision_detection_function == "SAT":
                if collides_optimized(polygon, polys[index2]):
                    collision_polygon = True
                    break
            elif collision_detection_function == "Naive":
                if collides(polygon, polys[index2]):
                    collision_polygon = True
                    break
        
        #if collides, fill it with shade. If no collision, don't shade it in 
        if save_plot or display_plot:
            if collision_polygon:
                ax.fill([vertex[0] for vertex in polygon], [vertex[1] for vertex in polygon], alpha=.25, fc='gray', ec='black')
            else:
                ax.fill([vertex[0] for vertex in polygon], [vertex[1] for vertex in polygon], alpha=.25, fc='white', ec='black')


    # Save the plot to a JPG file
    if save_plot:
        plt.savefig(output_file_name)
    
    if display_plot:
        plt.show()


#Algorithm to compute minkowski difference
#Given: Polygon P and Q
#Note: both polygons are an (n + 1) x 2 array. There are n vertices in the polygon but the 1st vertex has to be repeated at the end to indicate polygon is closed
#Optimized Version to Compute Minkowski Difference of P and Q
#CONFIRMED WORKS
def minkowski_difference_optimized(P, Q):
    #Goal is to compute Minowski Sum of P and -Q
    P, Q = P[:-1], Q[:-1]
    P, Q = P, -1 * Q 
        
    #Get starting pointers for P and Q
    P_y_values, Q_y_values = P[:, 1].flatten(), Q[:, 1].flatten()
    P_pointers, Q_pointers = np.argwhere(P_y_values == np.min(P_y_values)).flatten(), np.argwhere(Q_y_values == np.min(Q_y_values)).flatten()    
    P_pointer, Q_pointer = P_pointers[np.argmin(P[P_pointers, 0].flatten())], Q_pointers[np.argmin(Q[Q_pointers, 0].flatten())]
    
    #Initialize the counts for P and Q
    P_count, Q_count = 0,0
        
    #Initialize S
    S = []
    previous_P_phi, previous_Q_phi = None, None
                 
    #Iterate
    while P_count < P.shape[0] or Q_count < Q.shape[0]:
        #Sum up the two vertices and store it in S
        Pi, Qj = P[P_pointer % P.shape[0]], Q[Q_pointer % Q.shape[0]]   
        S.append(Pi + Qj)
            
        #Get edges Pi, P_i+1 and edge Qj, Q_j+1
        P_i1, Q_j1 = P[(P_pointer + 1) % P.shape[0]], Q[(Q_pointer + 1) % Q.shape[0]]
        P_edge, Q_edge = P_i1 - Pi, Q_j1 - Qj
            
        #Calculate polar angles. I don't want negative angles so just add 360 degrees. 
        P_phi, Q_phi = np.rad2deg(np.arctan2(P_edge[1], P_edge[0])), np.rad2deg(np.arctan2(Q_edge[1], Q_edge[0]))
        P_phi = P_phi + 360 if P_phi < 0 else P_phi
        Q_phi = Q_phi + 360 if Q_phi < 0 else Q_phi
            
        #Adjust if we have already made a full circle
        P_phi = P_phi + 360 if previous_P_phi != None and P_phi < previous_P_phi else P_phi
        Q_phi = Q_phi + 360 if previous_Q_phi != None and Q_phi < previous_Q_phi else Q_phi
                        
        #Adjust Pointer
        if P_phi < Q_phi:
            P_pointer, P_count = P_pointer + 1, P_count + 1
        elif P_phi > Q_phi:
            Q_pointer, Q_count = Q_pointer + 1, Q_count + 1
        else:
            P_pointer, P_count, Q_pointer, Q_count = P_pointer + 1, P_count + 1, Q_pointer + 1, Q_count + 1
            
        previous_P_phi, previous_Q_phi = P_phi, Q_phi
                         
    return np.array(S)

#Algorithm to compute minkowski difference
#Given: Polygon P and Q
#Note: both polygons are an (n + 1) x 2 array. There are n vertices in the polygon but the 1st vertex has to be repeated at the end to indicate polygon is closed
#CONFIRMED WORKS
def minkowski_difference(P, Q):
    P, Q = P[:-1], Q[:-1]
    S = []
        
    for vector1 in P:
        for vector2 in Q:
            S.append(vector1 - vector2)
        
    hull = ConvexHull(S)
    return hull.points[hull.vertices]


#Just a dummy test method. This method can be ignored :) 
#CONFIRMED WORKS
def compare_minkowskis(S, S2):
    # Use the sorted indices to rearrange the original array based on the first column
    sorted_indices = np.argsort(S[:, 0], axis = None)
    S = S[sorted_indices]
    
    # Use the sorted indices to rearrange the original array based on the first column
    sorted_indices = np.argsort(S2[:, 0], axis = None)
    S2 = S2[sorted_indices]

    return np.array_equal(S, S2)

#Formulating another Collision Detection Algorithm using Minkowski Differences
#Note: both polygons are an (n + 1) x 2 array. There are n vertices in the polygon but the 1st vertex has to be repeated at the end to indicate polygon is closed
#CONFIRMED WORKS
def collides_optimized_alternative(poly1, poly2):
    S = minkowski_difference_optimized(poly1, poly2)
    return point_contained(np.vstack((S, S[0])), np.array([0, 0]))