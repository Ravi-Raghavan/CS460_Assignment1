import numpy as np
import matplotlib.pyplot as plt
from collision_checking import*
from matplotlib.patches import Polygon, Circle, Rectangle
import math

# create arm (2 rectangular polygons with 2 joints, unlimited rotation) 
class PlanarArm:
    # constructor
    def __init__(self, ax):
        
        #initalize joint angles to 0 and radius of joint to 0.05
        self.theta1, self.theta2, self.radius = 0,0, 0.05

        # arm lengths and height:0.4m, 0.25m, 0.1m
        self.L1, self.L2, self.W =  0.4, 0.25, 0.1

        # load polygonal scene
        self.polygons = np.load('scene8.npy', allow_pickle= True)
        
        # set up the plot axis
        ax.set_xlim(0,2)
        ax.set_ylim(0,2)
        ax.set_aspect('equal')
        self.ax = ax
        self.fig = ax.figure

        # set up environment, remove any polygonals colliding with arm
        self.plot_arm()

        # removes polygons that collide with arm at start
       #polygon_collisions = self.collision()
       # if not polygon_collisions:
       #     for collision in polygon_collisions:
        #        self.polygons.remove(collision)

        # connect keyboard events to event handlers
        self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)

    # calculates positions and updates plot
    def plot_arm(self):
        # inital position
        x0,y0 = 1,1

        # calculate vertices, add 0.01 to the length to account for links (radius 0.05)
        x1 = x0 + (self.L1+.1) * np.cos(self.theta1) 
        y1 = x0 + (self.L1+.1) * np.sin(self.theta1)
        x2 = x1 + (self.L2+.1) * np.cos(self.theta2 + self.theta1) 
        y2 = y1 + (self.L2+.1) * np.sin(self.theta2 + self.theta1) 

        self.ax.clear()

        # generate arm joints
        self.generate_joint(self.radius, x0, y0)
        self.generate_joint(self.radius, x1, y1)
        self.generate_joint(self.radius, x2, y2)

        # generate arm links, second arm depends on angle of first
        box = np.array(self.generate_link((x0,y0),self.radius, self.theta1, self.L1, self.W))
        box2 = np.array(self.generate_link((x1,y1),self.radius, self.theta2+self.theta1, self.L2, self.W))

        # generate bounding boxes for joints (???)
        #square1 = Polygon(())

        # check for collisions
        self.collision(box)
        self.collision(box2)
        # if not empty, revert to previous move

        self.plot_polygons()
        self.ax.set_xlim(0,2)
        self.ax.set_ylim(0,2)

        #redraw
        self.fig.canvas.draw()

    # generate joints using radius and a given point
    def generate_joint(self, radius, x, y):
        circle = Circle((x,y),radius, fill = True, ec = 'blue', color = 'cornflowerblue')
        ax.add_patch(circle)

    # returns a rectangle ->represent rotating around a T
    def generate_link(self, joint, radius,  angle, length, width):
        # initialize rectangle points (increment by radius to account for joint)
        LL = (joint[0]+radius, joint[1]-width/2)
        LR = (joint[0]+radius+length, joint[1]-width/2)
        TR = (joint[0]+radius+length,joint[1]+width/2)
        TL = (joint[0]+radius,joint[1]+width/2)

        # get center point and shift points
        M = (joint[0],joint[1])
        LL = (LL[0]-M[0], LL[1]-M[1])
        LR = (LR[0]-M[0], LR[1]-M[1])
        TR =  (TR[0]-M[0], TR[1]-M[1])
        TL = (TL[0]-M[0], TL[1]-M[1])

        # apply rotation matrix
        LL = (LL[0]*np.cos(angle)-LL[1]*np.sin(angle), LL[0]*np.sin(angle)+LL[1]*np.cos(angle))
        LR = (LR[0]*np.cos(angle)-LR[1]*np.sin(angle), LR[0]*np.sin(angle)+LR[1]*np.cos(angle))
        TR = (TR[0]*np.cos(angle)-TR[1]*np.sin(angle), TR[0]*np.sin(angle)+TR[1]*np.cos(angle))
        TL = (TL[0]*np.cos(angle)-TL[1]*np.sin(angle), TL[0]*np.sin(angle)+TL[1]*np.cos(angle))

        # shift points back to get rectangle
        LL = (LL[0]+M[0], LL[1]+M[1])
        LR = (LR[0]+M[0], LR[1]+M[1])
        TR = (TR[0]+M[0], TR[1]+M[1])
        TL = (TL[0]+M[0], TL[1]+M[1])

        rectangle = [LL,LR,TR,TL, LL]
        arm1 = Polygon(rectangle, closed = True, ec = 'black', facecolor = 'cornflowerblue', alpha =0.5 )
        ax.add_patch(arm1)
        return rectangle

    # plots polygonal map
    def plot_polygons(self):
        for polygon in self.polygons:
             ax.add_patch(Polygon(polygon, closed = True, ec = 'black',facecolor = 'grey', alpha = 0.3))
        #redraw
        self.fig.canvas.draw()   

    # checks if bounding box collides with polygons
    def collision(self, box):
        obstacle = []
        for polygon in self.polygons:
            if collides(box, polygon):
                obstacle.append(polygon)
                print("Collision Detected: Previous position restored")
        return obstacle

    # responds press to key events, if a key causes a collision, it undos the step (how to access rectangles and circles???????)
    def on_key_press(self, event):
        # angle step
        step = np.pi/30
        if event.key == 'up':
            self.theta1 += step
        elif event.key == 'down':
            self.theta1 -= step
        elif event.key == 'right':
            self.theta2 += step
        elif event.key == 'left':
            self.theta2 -= step
        self.plot_arm()

# display plot
fig, ax = plt.subplots()
planarArm = PlanarArm(ax)
plt.show()