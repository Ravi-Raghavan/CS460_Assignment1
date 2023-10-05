from create_scene import *

# get user input 
P = int(input('number of polygons: '))
nmin = int(input('minimum vertices: '))
nmax = int(input('maximum vertices: '))
rmin = float(input('minimum radius: '))
rmax = float(input('maximum radius: ')) 
name = (input('name to save .npy file as: '))

constructScene(P,nmin,nmax,rmin,rmax,name)