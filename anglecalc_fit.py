import math
import numpy as np
from scipy.optimize import minimize

Qv_cal=np.array([1,1,1])
Qtheta_cal=np.array([3**(1/2),0,0])

        
# fitting process
def Omera_roration(omega):
    return np.array([[np.cos(np.radians(omega)),-np.sin(np.radians(omega)),0],[np.sin(np.radians(omega)),np.cos(np.radians(omega)),0],[0,0,1]])
def N_rotation(nu):
    return np.array([[1,0,0],[0,np.cos(np.radians(nu)),-np.sin(np.radians(nu))],[0,np.sin(np.radians(nu)),np.cos(np.radians(nu))]])
def M_rotation(mu):
    return np.array([[np.cos(np.radians(mu)),0,np.sin(np.radians(mu))],[0,1,0],[-np.sin(np.radians(mu)),0,np.cos(np.radians(mu))]])

def objective(angles):
    omega, mu, nu = angles
    rotation_matrix = Omera_roration(omega) @ N_rotation(mu) @ M_rotation(nu)
    transformed_vector = rotation_matrix @ Qv_cal
    return np.linalg.norm(transformed_vector - Qtheta_cal)

# 最適化によって omega, mu, nu を求める
initial_guess = [0, 0, 0]  # 初期値（全ての角度をゼロから開始）
result = minimize(objective, initial_guess, method='BFGS')

# 結果を表示
omega, mu, nu = result.x

print(f"omega: {omega}, mu: {mu}, nu: {nu}")