import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv

#open output file and make the list
file = open("./tests/output.csv", "r")
data = list(csv.reader(file, delimiter=","))
file.close()

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

#loop for update the image 
for k in range(0,len(data)-1,int(n)+1):
    #loop for create and update the bodies
    for j in range(int(n)+1):
        ax.scatter(float(data[k+j][1]),float(data[k+j][2]),float(data[k+j][3]))
    plt.draw()
    plt.pause(0.02)
    ax.cla()