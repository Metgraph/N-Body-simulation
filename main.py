import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv
import sys

#open output file and make the list
file = open(sys.argv[1], "r")
data = list(csv.reader(file, delimiter=","))
file.close()
scale_const=10e5

#detect how many bodies there are
n = float(0)
for i in range(len(data)):
    if n>float(data[i][0]):
        break
    n=float(data[i][0])

#create the 3d space


plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim3d(-scale_const, scale_const)
ax.set_ylim3d(-scale_const, scale_const)
ax.set_zlim3d(-scale_const, scale_const)


#loop for update the image 
for k in range(0,len(data)-1,int(n)+1):
    #loop for create and update the bodies
    for j in range(int(n)+1):
        ax.scatter(float(data[k+j][1]),float(data[k+j][2]),float(data[k+j][3]))
    plt.draw()
    plt.pause(0.02)
    ax.cla()
    ax.set_xlim3d(-scale_const, scale_const)
    ax.set_ylim3d(-scale_const, scale_const)
    ax.set_zlim3d(-scale_const, scale_const)
