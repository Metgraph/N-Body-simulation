import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv

file = open("../N-Body-simulation/tests/output.csv", "r")
data = list(csv.reader(file, delimiter=","))
file.close()

n = float(0)
for i in range(len(data)):
    if n>float(data[i][0]):
        break
    n=float(data[i][0])

plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for k in range(0,len(data)-1,int(n)+1):
    for j in range(int(n)+1):
        ax.scatter(float(data[k+j][1]),float(data[k+j][2]),float(data[k+j][3]))
    plt.draw()
    plt.pause(0.02)
    ax.cla()