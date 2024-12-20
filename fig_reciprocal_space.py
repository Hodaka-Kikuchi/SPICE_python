import configparser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d import Axes3D
import os
import sys

def plot_reciprocal_space(bpe, cphw, cp, fixe, sv1, sv2, **RLtable):
    """
    逆格子空間を描く（sp1vとsp2vがなす角度、格子点生成、ベクトルの表示）
    """
    astar = RLtable['astar']
    bstar = RLtable['bstar']
    cstar = RLtable['cstar']
    alpha_star = RLtable['alpha_star']
    beta_star = RLtable['beta_star']
    gamma_star = RLtable['gamma_star']
    n_a = RLtable['n_a']
    n_b = RLtable['n_b']
    n_c = RLtable['n_c']

    # sp1v, sp2vベクトルの計算
    sp1v = sv1[0] * astar + sv1[1] * bstar + sv1[2] * cstar
    sp2v = sv2[0] * astar + sv2[1] * bstar + sv2[2] * cstar
    
    # 角度計算 (cosθ = (sp1v ⋅ sp2v) / (||sp1v|| * ||sp2v||))
    dot_product = np.dot(sp1v, sp2v)
    norm_sp1v = np.linalg.norm(sp1v)
    norm_sp2v = np.linalg.norm(sp2v)
    
    cosine_theta = dot_product / (norm_sp1v * norm_sp2v)
    angle_deg = np.degrees(np.arccos(np.clip(cosine_theta, -1.0, 1.0)))  # 角度を度で計算

    print(f"Angle between sp1v and sp2v: {angle_deg:.2f} degrees")

    # 格子点の計算 (整数倍で描画)
    n_points = 10  # 格子点の数
    grid_points = []

    # 格子点の生成
    for i in range(-n_points, n_points+1):
        for j in range(-n_points, n_points+1):
            grid_point = i * sp1v + j * sp2v
            grid_points.append(grid_point)

    grid_points = np.array(grid_points)
    
    # hkl_calベクトルの計算
    hkl_cal = cp[0] * astar + cp[1] * bstar + cp[2] * cstar
    
    # 2D空間にプロット
    sp1v_2d = sp1v[:2]  # x, y成分のみ
    sp2v_2d = sp2v[:2]  # x, y成分のみ
    grid_points_2d = grid_points[:, :2]  # x, y成分のみ
    hkl_cal_2d = hkl_cal[:2]  # hkl_calもx, y成分のみ

    # プロットの準備
    plt.figure(figsize=(6, 6))
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    
    # 格子点を描画
    plt.scatter(grid_points_2d[:, 0], grid_points_2d[:, 1], color='blue', s=10, label='Lattice Points')
    
    # hkl_calベクトルを描画
    plt.quiver(0, 0, hkl_cal_2d[0], hkl_cal_2d[1], angles='xy', scale_units='xy', scale=1, color='red', label='hkl_cal')
    
    # sp1v, sp2vベクトルを描画
    plt.quiver(0, 0, sp1v_2d[0], sp1v_2d[1], angles='xy', scale_units='xy', scale=1, color='green', label='sp1v')
    plt.quiver(0, 0, sp2v_2d[0], sp2v_2d[1], angles='xy', scale_units='xy', scale=1, color='orange', label='sp2v')
    
    # 軸の設定
    plt.xlim(-1.5 * max(norm_sp1v, norm_sp2v), 1.5 * max(norm_sp1v, norm_sp2v))
    plt.ylim(-1.5 * max(norm_sp1v, norm_sp2v), 1.5 * max(norm_sp1v, norm_sp2v))
    
    plt.gca().set_aspect('equal', adjustable='box')
    
    # ラベルとタイトル
    plt.title('Reciprocal Space Plot (2D)')
    plt.xlabel('k1')
    plt.ylabel('k2')
    
    plt.legend()

    #plt.show()
