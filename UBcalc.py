# UBマトリックスの計算。
import math
import numpy as np

def UB_calc(sp1,sp2,astar,bstar,cstar):
    u1 = sp1[0]*astar+sp1[1]*bstar+sp1[2]*cstar
    U1 = u1 / np.sqrt(np.dot(u1, u1))
    u2 = sp2[0]*astar+sp2[1]*bstar+sp2[2]*cstar
    uu2 = u2 - np.dot(U1, u2) * U1
    U2 = uu2 / np.linalg.norm(uu2)
    U3=[U1[1]*U2[2]-U1[2]*U2[1], U1[2]*U2[0]-U1[0]*U2[2], U1[0]*U2[1]-U1[1]*U2[0]]
    
    # 結果を返す
    return U1,U2,U3
    