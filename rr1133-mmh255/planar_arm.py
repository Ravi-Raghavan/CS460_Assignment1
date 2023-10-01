import numpy as np
import matplotlib.pyplot as plt
from collision_checking import*
from matplotlib.colors import ListedColormap
from matplotlib.patches import Polygon, Circle
import argparse

class PlanarArm:
    def __init__(self, ax, mode):
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
        self.ax = ax
        self.fig = ax.figure
        ax.set_xlim(0,2)
        ax.set_ylim(0,2)
        ax.set_aspect('equal')

        if mode == True:
            # set up environment with initial arm at (1,1) 
            self.plot_arm()
        else:
            # display cspace for the arm and polygonal map
            self.occupancy_grid()

        # connect keyboard events to event handlers
        self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)

    # plots the 2R arm alongside polygonal map
    def plot_arm(self):
        self.ax.clear()
        self.robot.clear() # reset robot of previous position and orientation 

        # store vertices of robot
        points = self.generate_robot()

        # generate arm joints using vertices & adds to plot
        self.generate_joint(self.radius, points[0], points[1])
        self.generate_joint(self.radius, points[2], points[3])
        self.generate_joint(self.radius, points[4], points[5])

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

    def generate_robot(self):
        # inital position
        x0,y0 = 1,1

        # calculate vertices, add 0.1 to the length to account for links (radius 0.05)
        x1 = x0 + (self.L1+.1) * np.cos(self.theta1) 
        y1 = y0 + (self.L1+.1) * np.sin(self.theta1)
        x2 = x1 + (self.L2+.1) * np.cos(self.theta2 + self.theta1) 
        y2 = y1 + (self.L2+.1) * np.sin(self.theta2 + self.theta1) 

        # construct parts of robot
        self.robot.append(np.array(self.get_link((x0,y0),self.radius, self.theta1, self.L1, self.W)))
        self.robot.append(np.array(self.get_link((x1,y1),self.radius, self.theta2+self.theta1, self.L2, self.W)))
        self.robot.append(np.array(self.get_link((x0,y0),-0.05, self.theta1, self.radius*2, self.radius*2)))
        self.robot.append(np.array(self.get_link((x1,y1), -0.05, self.theta2, self.radius*2, self.radius*2)))
        self.robot.append(np.array(self.get_link((x2,y2), -0.05, self.theta2, self.radius*2, self.radius*2)))

        return [x0,y0, x1, y1, x2, y2]

    # generate a joint and add it to the plot
    def generate_joint(self, radius, x, y):
        circle = Circle((x,y),radius, fill = True, ec = 'blue', color = 'cornflowerblue')
        self.ax.add_patch(circle)

    # generate a link (doesn't add to plot if it is a box for a joint)
    def generate_link(self, rectangle, boundBox):
        if boundBox:
            arm1 = Polygon(rectangle)
        else:
            arm1 = Polygon(rectangle, closed = True, ec = 'black', facecolor = 'cornflowerblue', alpha =0.5 )
            self.ax.add_patch(arm1)

    # creates a link
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
        for polygon in self.polygons:
            for part in self.robot:
                if collides_optimized(part, polygon):
                    collision_detected = True
                    break
        
        # reverts angles and generates previous arm
        if collision_detected:
            self.theta1 = prev_theta1
            self.theta2 = prev_theta2
            self.plot_arm() 

    # increments/decrements theta1 and theta2 depending on key 
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

        # check for collisions (revert if one occurs)
        self.check_collisions(prev_theta1, prev_theta2)

    # compute and visualize cspace of workspace
    def occupancy_grid(self):
        data = np.zeros((100,100), dtype = int)
        step = np.pi/50

        # for every increment of theta1, theta2 does a full rotation
        for x in range(100):
            self.theta2 = 0.0
            for y in range(100):
                detected = False
                # generate a new robot after each increment
                self.robot.clear()
                self.generate_robot()
                # check for collisions
                for polygon in self.polygons:
                    for part in self.robot:
                        if collides_optimized(part,polygon):
                            data[x][y] = 1
                            # if collision occured, break to check collision with other polygons
                            detected = True
                            if detected:
                                 break
                self.theta2 += step
            self.theta1 += step

        # visualize 
        cspace = np.transpose(data) #adjust 2D array so the plot orientation is correct
        plt.figure(figsize=(7, 7))
        plt.imshow(cspace, extent= (0,100, 100,0), cmap = ListedColormap(['mediumpurple', 'yellow']))
        plt.title('Collision-Free Configuration Space')
        plt.xlabel('Joint Angle 1')
        plt.ylabel('Joint Angle 2')
        plt.show()

def main():
     # create command-line argument parser
    parser = argparse.ArgumentParser(description = 'Planar Arm Program')
    
    # add arguments
    parser.add_argument('--mode', choices=['arm', 'cspace'], default='arm',
                        help='Specify the initial mode (arm or cspace)')
    
    # Parse the command-line arguments
    args = parser.parse_args()

    # display and control arm OR display configuration space depending on argument
    fig, ax = plt.subplots()
    if args.mode == 'arm':
        robot = PlanarArm(ax, True)
    elif args.mode == 'cspace':
        robot = PlanarArm(ax, False)
    plt.show()

if __name__ == "__main__":
    main()