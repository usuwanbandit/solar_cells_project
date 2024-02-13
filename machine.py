import numpy as np
from numpy import *
import cmath
Z = 48+64j
I = 2.5

V = I*Z

P = V* np.conj(I)
print("Z", Z, "abs(Z)", abs(Z))
print("I", I, "abs(I)", abs(I))
print("V", V ,"abs(V)", abs(V))
print("P", P, "abs(P)", abs(P))
print(P)
print(P.real)
# i = np.sqrt(-1)
# Z = 625 + 75*i
# print(Z)