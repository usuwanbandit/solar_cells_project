import numpy as np
from numpy import *
import matplotlib.pyplot as plt


def find_slope(list1, list2):
    # Check if both lists have the same length
    if len(list1) != len(list2):
        raise ValueError("Lists must have the same length.")

    # Calculate the differences between corresponding elements
    differences_x = [list1[i + 1] - list1[i] for i in range(len(list1) - 1)]
    differences_y = [list2[i + 1] - list2[i] for i in range(len(list2) - 1)]

    # Calculate the slope
    slope = differences_y[-1] / differences_x[-1]

    return slope


V1 = [220, 200, 160, 120, 80, 40, 0]
V2 = [121.8, 111.3, 89, 66.4, 43.52, 20.1, 0]

A1 = [0.9, 0.75, 0.6, 0.45, 0.3, 0.15, 0]
A2 = [1.65, 1.4, 1.2, 0.76, 0.52, 0.24, 0]

V2_4 = [121, 120, 119, 118, 116, 115, 110]
A2_4 = [0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8]
plt.plot(A2_4, V2_4)
plt.xlabel("A(A)")
plt.ylabel("V(V)")
plt.title(f" VR = (V$_{'NL'}$ - V$_{'FL'}$)/V$_{'FL'}$ * 100 = {(121.8-110)/110 *100} ")
plt.legend()



plt.show()
