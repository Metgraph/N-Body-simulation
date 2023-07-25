from numpy._typing import NDArray
import matplotlib.pyplot as plt
import numpy as np
import csv
import sys

def open_file(file_name: str) -> tuple[NDArray, int]:
    print(f"Opening file: {file_name}, it might take a while.", end="\r")

    with open(file_name, "r") as f:
        data = list(csv.reader(f, delimiter=","))

    data = np.array([list(map(float, x)) for x in data])

    print("Read completed!                                          ")
    return data, int(np.max(data[:, :1])) + 1


data1, n_corps1 = open_file(sys.argv[1])
data2, n_corps2 = open_file(sys.argv[2])

if n_corps1 != n_corps2:
    print("Il numero dei corpi non corrisponde")
    exit(1)

diff = np.absolute(data1 - data2)

i = np.arange(len(diff)) % n_corps1 == 0



y = diff[i]
y = y[30:101]
x = np.arange(len(y))


coefficients = np.polyfit(x, y[:, 1], 1)
m, b = coefficients

# Calcola i valori della retta di regressione
regression_line = m * x + b
#fitted_curve = np.polyval(coefficients, x)

print(y)

plt.plot(x, y[:, 1:2], 'bo-', linewidth=0.5, markersize=2)
plt.plot(x, regression_line, 'r-', linewidth=1)  # Aggiungi la retta di regressione rossa
plt.xscale('linear')
#plt.yscale('log')
plt.yscale('linear')
plt.xlabel('t')
plt.ylabel('Diff x')
plt.title('Plot in scala logaritmica')
plt.grid(True)
plt.show()


