import numpy as np
import matplotlib.pyplot as plt

# --- 楕円係数の導出関数 ---
def ellipse_coefficients(RM, plane=("Qx", "E")):
    # 軸のインデックス
    axis_map = {"Qx": 0, "Qy": 1, "E": 2, "Qz": 3}
    i = axis_map[plane[0]]
    j = axis_map[plane[1]]

    # 残す軸 i, j の係数
    A = RM[i, i]
    C = RM[j, j]
    B = 2 * RM[i, j]

    # 残りの軸（消すべき軸）を探す
    other_axes = [k for k in range(4) if k not in [i, j]]

    # ガウス分布の周辺化に基づく修正項（Schur補行列的なもの）
    for k in other_axes:
        A -= (RM[i, k] ** 2) / RM[k, k]
        C -= (RM[j, k] ** 2) / RM[k, k]
        B -= 2 * RM[i, k] * RM[j, k] / RM[k, k]

    # 楕円の式の右辺定数項
    F = -2 * np.log(2)
    
    return A, B, C, F

# --- 楕円のプロット関数 ---
def plot_ellipse_from_coeff(A, B, C, F, xlim, ylim, xlabel, ylabel, color):
    x = np.linspace(-xlim, xlim, 500)
    y = np.linspace(-ylim, ylim, 500)
    X, Y = np.meshgrid(x, y)
    
    ellipse = A * X**2 + B * X * Y + C * Y**2 + F

    plt.contour(X, Y, ellipse, levels=[0], colors=color)
    plt.axhline(0, color='black', linestyle='--', linewidth=0.5)
    plt.axvline(0, color='black', linestyle='--', linewidth=0.5)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.axis('equal')

# --- 共分散行列（例） ---
RM = np.array([
    [48286,   3178,    -2992,     -0.0],
    [3178,  21203, -3175,     0.0],
    [-2992, -3175,  710,     -0.0],
    [-0.0,     0.0,     -0.0,    1551]
])

# --- Qx vs E 楕円 ---
A1, B1, C1, F1 = ellipse_coefficients(RM, plane=("Qx", "E"))
plt.figure(figsize=(6,5))
plot_ellipse_from_coeff(A1, B1, C1, F1, xlim=0.05, ylim=1.0, xlabel="Qx (Å⁻¹)", ylabel="E (meV)", color="red")
plt.title("Qx vs E Resolution Ellipse")
plt.show()

# --- Qy vs E 楕円 ---
A2, B2, C2, F2 = ellipse_coefficients(RM, plane=("Qy", "E"))
plt.figure(figsize=(6,5))
plot_ellipse_from_coeff(A2, B2, C2, F2, xlim=0.05, ylim=1.0, xlabel="Qy (Å⁻¹)", ylabel="E (meV)", color="blue")
plt.title("Qy vs E Resolution Ellipse")
plt.show()
