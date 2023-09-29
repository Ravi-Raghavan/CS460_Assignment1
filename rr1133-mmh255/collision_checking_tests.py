import numpy as np 
import matplotlib.pyplot as plt
from collision_checking import *

#Testing this on multiple polygon scenes
scenes = ["scene.npy", "scene2.npy", "scene3.npy", "scene4.npy", "scene5.npy", "scene6.npy", "scene7.npy", "scene8.npy"]
scenes = ["scene.npy", "scene2.npy", "scene3.npy", "scene4.npy", "scene5.npy", "scene6.npy", "scene7.npy"]
for index, scene in enumerate(scenes):
    if index < 6: continue
    polygons = np.load(f"rr1133-mmh255/{scene}", allow_pickle = True)
    file_name = f"rr1133-mmh255/Problem2_scene{index + 1}.jpg"
    print_diagnostics = False
    if (index == 1):
        print_diagnostics = True
        
    plot(polygons, file_name, False, print_diagnostics)