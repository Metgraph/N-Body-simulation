import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv
import sys

"""
Create first a python environment: $> python3 -m venv venv
Activate it: $> source venv/bin/activate
Install dependencies: $> python3 -m pip install -r requirements.txt

Exit from environment: $> deactivate
Update dependencies (overwrite requirements file): $> pip freeze > requirements.txt
"""

scale=10e10
pause_time = 0.02

def open_file():
    with open(sys.argv[1], "r") as f:
        data = list(csv.reader(f, delimiter=","))

    # convert data and check numbers of bodies
    n = 0
    for i in range(len(data)):
        for j in range(len(data[i])):
            if j == 0:
                data[i][j] = int(data[i][j])
                if data[i][j] > n:
                    n = data[i][j]
            else:
                data[i][j] = float(data[i][j])
    return data, n

def show_plot(data, bodies, scale, pause_time):
    #create the 3d space
    plt.ion()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim3d(-scale, scale)
    ax.set_ylim3d(-scale, scale)
    ax.set_zlim3d(-scale, scale)

    #loop for update the image
    for k in range(0,len(data)-1, bodies+1):
        #loop for create and update the bodies
        for j in range(bodies + 1):
            ax.scatter(data[k+j][1], data[k+j][2], data[k+j][3])
        plt.draw()
        plt.pause(pause_time)
        ax.cla()
        ax.set_xlim3d(-scale, scale)
        ax.set_ylim3d(-scale, scale)
        ax.set_zlim3d(-scale, scale)

def main():
    global scale
    global pause_time

    argc = len(sys.argv)
    if argc < 2:
        print(f'Usage: python3 {sys.argv[0]} FILE [SCALE|<default> [PAUSE_TIME]]')
        print(f'Default scale: {scale}')
        print(f'Default pause time: {pause_time}')
        exit(1)

    i = 2
    while i < argc:
        if i == 2:
            if sys.argv[2][0].isnumeric():
                scale = float(sys.argv[2])
        elif i == 3:
            pause_time = float(sys.argv[3])
        i += 1

    print(f'Using scale: {scale}')
    print(f'Using pause time: {pause_time}')

    data, bodies = open_file()
    show_plot(data, bodies, scale, pause_time)

if __name__ == "__main__":
    main()
