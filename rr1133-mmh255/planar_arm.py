import numpy as np
import matplotlib.pyplot as plt
from collision_checking import*
from matplotlib.patches import Polygon, Circle

class PlanarArm:
    def __init__(self, ax):
        # initalize joint angles 
        self.theta1, self.theta2, self.radius = 0,0, 0.05

        # arm lengths and height
        self.L1, self.L2, self.W =  0.4, 0.25, 0.1

        # keeps track of all parts of robot
        self.robot = []

        # load polygonal scene
        self.polygons = np.load('arm_polygons.npy', allow_pickle = True)
        self.start = True
        
        # set up the plot axis
        ax.set_xlim(0,2)
        ax.set_ylim(0,2)
        ax.set_aspect('equal')
        self.ax = ax
        self.fig = ax.figure

        # set up environment with initial arm at (1,1) at 0 degrees
        self.plot_arm()

        # connect keyboard events to event handlers
        self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)

    # calculates position and orientation of arm
    def plot_arm(self):
        self.ax.clear()
        self.robot.clear()

        # inital position
        x0,y0 = 1,1

        # calculate vertices, add 0.1 to the length to account for links (radius 0.05)
        x1 = x0 + (self.L1+.1) * np.cos(self.theta1) 
        y1 = y0 + (self.L1+.1) * np.sin(self.theta1)
        x2 = x1 + (self.L2+.1) * np.cos(self.theta2 + self.theta1) 
        y2 = y1 + (self.L2+.1) * np.sin(self.theta2 + self.theta1) 

        # generate arm joints
        self.generate_joint(self.radius, x0, y0)
        self.generate_joint(self.radius, x1, y1)
        self.generate_joint(self.radius, x2, y2)

        # generate arm links
        self.robot.append(np.array(self.get_link((x0,y0),self.radius, self.theta1, self.L1, self.W)))
        self.robot.append(np.array(self.get_link((x1,y1),self.radius, self.theta2+self.theta1, self.L2, self.W)))

        # generate bounding boxes for joints 
        self.robot.append(np.array(self.get_link((x0,y0),-0.05, self.theta1, self.radius*2, self.radius*2)))
        self.robot.append(np.array(self.get_link((x1,y1), -0.05, self.theta2, self.radius*2, self.radius*2)))
        self.robot.append(np.array(self.get_link((x2,y2), -0.05, self.theta2, self.radius*2, self.radius*2)))

        # adds links to plot
        self.generate_link(np.array(self.robot[0]), False)
        self.generate_link(np.array(self.robot[1]), False)
        self.generate_link(np.array(self.robot[2]), True)
        self.generate_link(np.array(self.robot[3]), True)
        self.generate_link(np.array(self.robot[4]), True)

        # draw 
        self.plot_polygons()
        self.ax.set_xlim(0,2)
        self.ax.set_ylim(0,2)
        self.fig.canvas.draw()

    # plots polygonal map
    def plot_polygons(self):
        # removes any polygons colliding with arm at START ENVIRONMENT 
        if self.start:
            for part in self.robot:
                self.start_collision(part)   
        self.start = False

        # adds polygons to map
        for polygon in self.polygons:
             self.ax.add_patch(Polygon(polygon, closed = True, ec = 'black',facecolor = 'grey', alpha = 0.3))

        #redraw
        self.fig.canvas.draw()  

    # generate a joint and add it to the plot
    def generate_joint(self, radius, x, y):
        circle = Circle((x,y),radius, fill = True, ec = 'blue', color = 'cornflowerblue')
        self.ax.add_patch(circle)

    # generate a link (doesn't add to plot if it is a joint box)
    def generate_link(self, rectangle, boundBox):
        if boundBox:
            arm1 = Polygon(rectangle)
        else:
            arm1 = Polygon(rectangle, closed = True, ec = 'black', facecolor = 'cornflowerblue', alpha =0.5 )
            self.ax.add_patch(arm1)

    # returns a rectangle -> represent rotating around a T
    def get_link(self, joint, radius,  angle, length, width):
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

        return rectangle 

    # checks for any start collisions and removes them
    def start_collision(self, part):
        to_delete = []
        for i, polygon in enumerate(self.polygons):
            if collides_optimized(part, polygon):
                to_delete.append(i)
        self.polygons= np.delete(self.polygons, to_delete, 0)

    # checks if polygon collision is detected (if yes, revert robot)
    def check_collisions(self, prev_theta1, prev_theta2):
        # constructs new arm after key event
        self.plot_arm()

        # checks if any part of arm collides with polygon
        collision_detected = False
        for part in self.robot:
            for polygon in self.polygons:
                if collides_optimized(part, polygon):
                    collision_detected = True
                    break
        
        # reverts angles and generates previous arm
        if collision_detected:
            self.theta1 = prev_theta1
            self.theta2 = prev_theta2
            self.plot_arm() 

    def on_key_press(self, event):
        # hold previous joint angles
        prev_theta1 = self.theta1
        prev_theta2 = self.theta2

        step = np.pi/30 #angle step
        self.event = event.key
        if event.key == 'up':
            self.theta1 += step
        elif event.key == 'down':
            self.theta1 -= step
        elif event.key == 'right':
            self.theta2 += step
        elif event.key == 'left':
            self.theta2 -= step

        # check for collisions 
        self.check_collisions(prev_theta1, prev_theta2)

    # use start environment to visualize the workspace? configuration space?
   # def occupancy_grid(self):
        # divide 2pi/100
        # 100 x 100 grid, each cell = value of 2 joint angles
        # YELLOW = collisions with polygonal obstacles
        # x = first joint
        # y = second joint
        # how to create grid filled with color??? 
        # get


# display plot
fig, ax = plt.subplots()
planarArm = PlanarArm(ax)
plt.show()