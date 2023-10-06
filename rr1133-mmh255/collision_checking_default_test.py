import numpy as np 
import matplotlib.pyplot as plt
from collision_checking import *

polygons = np.load(f"collision_checking_polygons.npy", allow_pickle = True)
file_name = f"Problem2_defaultScene.jpg"

plot(polygons, file_name, False, True, False, collision_detection_function = "Minkowski")