import configparser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d import Axes3D
import os
import sys

def solve_linear_combination(sp1v, sp2v, hkl_cal):
    """
    hkl_cal = c1 * sp1v + c2 * sp2v の線形方程式を解いて c1, c2 を求める
    """
    # sp1v, sp2v を列ベクトルとして並べた行列を作成
    A = np.vstack([sp1v, sp2v]).T  # 2行3列の行列

    # hkl_cal を列ベクトルとして用意
    b = hkl_cal

    # A * [c1, c2] = b の方程式を解く
    coefficients = np.linalg.lstsq(A, b, rcond=None)[0]  # 最小二乗法で解く

    c1, c2 = coefficients
    return c1, c2

def plot_reciprocal_space(bpe, bpc2, cphw, cp, fixe, sv1, sv2, RLtable, angletable):
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
    
    if fixe==0: # ei fix
        Ei = bpe
        Ef = bpe - cphw
    elif fixe==1: # ef fix
        Ei = bpe + cphw
        Ef = bpe
    ki_cal=(Ei/2.072)**(1/2)
    kf_cal=(Ef/2.072)**(1/2)
    
    # sp1v, sp2vベクトルの計算
    sp1v = sv1[0] * astar + sv1[1] * bstar + sv1[2] * cstar
    sp2v = sv2[0] * astar + sv2[1] * bstar + sv2[2] * cstar
    
    # 角度計算 (cosθ = (sp1v ⋅ sp2v) / (||sp1v|| * ||sp2v||))
    dot_product = np.dot(sp1v, sp2v)
    norm_sp1v = np.linalg.norm(sp1v)
    norm_sp2v = np.linalg.norm(sp2v)
    
    cosine_theta = dot_product / (norm_sp1v * norm_sp2v)
    angle_deg = np.degrees(np.arccos(np.clip(cosine_theta, -1.0, 1.0)))  # 角度を度で計算

    # XY空間上にベクトルを投影
    sp1v_2d = [norm_sp1v,0]
    sp2v_2d = [norm_sp2v * np.cos(np.radians(angle_deg)),norm_sp2v * np.sin(np.radians(angle_deg))]

    # 格子点の計算 (整数倍で描画)
    n_points1 = int(ki_cal/norm_sp1v)+1  # 格子点の数
    grid_points = []

     # 格子点の生成 (i, j は整数)
    n_points2 = int(ki_cal/norm_sp2v)+1  # 格子点の数
    grid_points = []

    # 格子点の生成 (sp1v_2d と sp2v_2d の線形和)
    for i in range(-n_points1, n_points1+1):
        for j in range(-n_points2, n_points2+1):
            grid_point = i * np.array(sp1v_2d) + j * np.array(sp2v_2d)  # 格子点は2次元ベクトルとして計算
            grid_points.append(grid_point)

    grid_points = np.array(grid_points)  # 格子点をnumpy配列に変換
    
    # hkl_calベクトルの計算 (sp1v_2d と sp2v_2d の線形和)
    hkl_cal = cp[0] * astar + cp[1] * bstar + cp[2] * cstar
    
    # 線形結合を解く
    c1, c2 = solve_linear_combination(sp1v, sp2v, hkl_cal)
    
    # プロットの準備
    plt.figure(figsize=(8, 5))
    plt.subplots_adjust(left=0.01, bottom=0.25)
    
    # 格子点を描画
    plt.scatter(grid_points[:, 0], grid_points[:, 1], color='black', s=10, label='Lattice Points')
    
    hkl_x = c1 * sp1v_2d[0] + c2 * sp2v_2d[0]
    hkl_y = c1 * sp1v_2d[1] + c2 * sp2v_2d[1]
    # hkl_calベクトルを描画
    plt.quiver(0, 0, hkl_x, hkl_y, angles='xy', scale_units='xy', scale=1, color='red', label='hkl_cal')
    
    # sp1v_2d, sp2v_2dベクトルを描画
    # 端に表示するようにする。
    plt.quiver(-sp1v_2d[0]*n_points1-sp2v_2d[0]*n_points2, -sp1v_2d[1]*n_points1-sp2v_2d[1]*n_points2, sp1v_2d[0], sp1v_2d[1], angles='xy', scale_units='xy', scale=1, color='magenta', label='axis 1')
    plt.quiver(-sp1v_2d[0]*n_points1-sp2v_2d[0]*n_points2, -sp1v_2d[1]*n_points1-sp2v_2d[1]*n_points2, sp2v_2d[0], sp2v_2d[1], angles='xy', scale_units='xy', scale=1, color='aqua', label='axis 2')
    
    # ki, kfベクトルを描画. ki//y_axis
    #inst_x = ki_cal * np.sin(np.radians(angletable['C2']-angletable['offset']))
    #inst_y = ki_cal * np.cos(np.radians(angletable['C2']-angletable['offset']))
    ki_vx = ki_cal * np.sin(np.radians(angletable['C2']-angletable['offset']))
    ki_vy = ki_cal * np.cos(np.radians(angletable['C2']-angletable['offset']))
    
    kf_vx = kf_cal * np.sin(np.radians(angletable['C2']-angletable['offset']-angletable['A2']))
    kf_vy = kf_cal * np.cos(np.radians(angletable['C2']-angletable['offset']-angletable['A2']))
    plt.quiver(hkl_x-ki_vx, hkl_y-ki_vy, ki_vx, ki_vy,  angles='xy', scale_units='xy', scale=1, color='blue', label='ki')
    plt.quiver(hkl_x-ki_vx, hkl_y-ki_vy, kf_vx, kf_vy, angles='xy', scale_units='xy', scale=1, color='green', label='kf')
    
    # 軸の設定
    #plt.xlim(-1.5 * max(norm_sp1v, norm_sp2v), 1.5 * max(norm_sp1v, norm_sp2v))
    #plt.ylim(-1.5 * max(norm_sp1v, norm_sp2v), 1.5 * max(norm_sp1v, norm_sp2v))
    
    plt.gca().set_aspect('equal', adjustable='box')
    
    # ラベルとタイトル
    plt.title('Reciprocal Space Plot')
    plt.xlabel(r'$k_x\ (\mathrm{\AA}^{-1})$')
    plt.ylabel(r'$k_y\ (\mathrm{\AA}^{-1})$')
    
    # 凡例をグラフの外に配置
    plt.legend(loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0., fontsize=10)

    #plt.show()
