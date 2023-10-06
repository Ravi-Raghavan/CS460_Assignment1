import numpy as np 
import matplotlib.pyplot as plt
from collision_checking import *
import time

#Testing this on multiple polygon scenes
scenes = ["scene1.npy", "scene2.npy", "scene3.npy", "scene4.npy", "scene5.npy", "scene6.npy", "scene7.npy", "scene8.npy"]

minkowski_array = []
sat_array = []
naive_array = []

for index, scene in enumerate(scenes):
    polygons = np.load(f"p2_{scene}", allow_pickle = True)
    file_name = f"Problem2_scene{index + 1}.jpg"
    print_diagnostics = False
    
    minkowski_time = 0
    sat_time = 0
    naive_time = 0
            
    #Run 10 Trials to Average Results and get some average value
    cumulative_time = 0
    trials = 10
    for trial in range(trials):
        start_time = time.time()
        plot(polygons, file_name, False, False, print_diagnostics, collision_detection_function = "Minkowski")
        end_time = time.time()
        cumulative_time += (end_time - start_time)
    
    minkowski_time = (cumulative_time) / trials
    
    #Run 10 Trials to Average Results and get some average value
    cumulative_time = 0
    trials = 10
    for trial in range(trials):
        start_time = time.time()
        plot(polygons, file_name, False, False, print_diagnostics, collision_detection_function = "SAT")
        end_time = time.time()
        cumulative_time += (end_time - start_time)
    
    sat_time = (cumulative_time) / trials
    
    #Run 10 Trials to Average Results and get some average value
    cumulative_time = 0
    trials = 10
    for trial in range(trials):
        start_time = time.time()
        plot(polygons, file_name, False, False, print_diagnostics, collision_detection_function = "Naive")
        end_time = time.time()
        cumulative_time += (end_time - start_time)
    
    naive_time = (cumulative_time) / trials
    
    # print(f"Scene #: {index + 1}, Minkowski Time: {minkowski_time:.10f}, SAT Time: {sat_time:.10f}, Naive Time: {naive_time:.10f}")
    print(f"{index + 1} & {naive_time:.10f} & {sat_time} & {minkowski_time:.10f} \\\\")
    
    #Append to arrays
    minkowski_array.append(minkowski_time)
    sat_array.append(sat_time)
    naive_array.append(naive_time)

minkowski_vs_naive = 0
sat_vs_naive = 0
sat_vs_minkowski = 0

for index in range(len(scenes)):
    if minkowski_array[index] < naive_array[index]:
        minkowski_vs_naive += 1
        
    if sat_array[index] < naive_array[index]:
        sat_vs_naive += 1
    
    if sat_array[index] < minkowski_array[index]:
        sat_vs_minkowski += 1


print(f"Minkowski Out Performs Naive in {minkowski_vs_naive} trials")   
print(f"SAT Out Performs Naive in {sat_vs_naive} trials")
print(f"SAT Out Performs Minkowski in {sat_vs_minkowski} trials")   