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

def plot_reciprocal_space(bpe, bpc2, cphw, cp, fixe, sv1, sv2, RLtable,A_sets, C_sets,QE_sets, initial_index=0):
    """
    逆格子空間を描く（sp1vとsp2vがなす角度、格子点生成、ベクトルの表示）
    """
    #QE_sets.append([hw, h, k, l])
    
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
        Ef = bpe - QE_sets[initial_index][0]
    elif fixe==1: # ef fix
        Ei = bpe + QE_sets[initial_index][0]
        Ef = bpe
    ki_cal=(Ei/2.072)**(1/2)
    kf_cal=(Ef/2.072)**(1/2)
    
    """
    # ベクトル u1, u2 の計算
    u1 = sv1[0] * astar + sv1[1] * bstar + sv1[2] * cstar
    U1 = u1 / np.linalg.norm(u1)
    
    u2 = sv2[0] * astar + sv2[1] * bstar + sv2[2] * cstar
    uu2 = u2 - np.dot(U1, u2) * U1
    U2 = uu2 / np.linalg.norm(uu2)
    """
    
    # --- 3D格子ベクトルを計算 ---
    sp1v = sv1[0] * astar + sv1[1] * bstar + sv1[2] * cstar
    sp2v = sv2[0] * astar + sv2[1] * bstar + sv2[2] * cstar

    # --- 2D基底を構築 ---
    # e1: sp1v を正規化したベクトル (x軸方向)
    e1 = sp1v / np.linalg.norm(sp1v)

    # e2: sp2v を e1 に直交化 (Gram-Schmidt法)
    proj = np.dot(sp2v, e1) * e1
    orth = sp2v - proj
    e2 = orth / np.linalg.norm(orth)

    # --- 3D→2D変換行列 ---
    basis2D = np.vstack([e1, e2])  # shape (2,3)

    def project_to_2d(vec3):
        """3Dベクトルを2D平面に射影"""
        return basis2D @ vec3

    # --- 2D座標に変換 ---
    sp1v_2d = project_to_2d(sp1v)
    sp2v_2d = project_to_2d(sp2v)

    # --- 格子点生成 ---
    n_points1 = int(2*ki_cal/np.linalg.norm(sp1v))+1
    n_points2 = int(2*ki_cal/np.linalg.norm(sp2v))+1

    grid_points = []
    for i in range(-n_points1, n_points1+1):
        for j in range(-n_points2, n_points2+1):
            grid_points.append(i*sp1v_2d + j*sp2v_2d)

    grid_points = np.array(grid_points)

    # --- hklベクトルを2Dへ ---
    hkl_cal = QE_sets[initial_index][1]*astar + QE_sets[initial_index][2]*bstar + QE_sets[initial_index][3]*cstar
    hkl_cal_2d = project_to_2d(hkl_cal)

    # --- プロット ---
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(grid_points[:,0], grid_points[:,1], s=10, color="blue", alpha=0.5)
    ax.scatter(hkl_cal_2d[0], hkl_cal_2d[1], color="red", marker="x", s=100)
    
    # スライダー設定
    ax_slider = plt.axes([0.25, 0.10, 0.65, 0.03], facecolor='lightgoldenrodyellow')
    slider = Slider(ax_slider, 'scan number', 1, len(QE_sets) , valinit=initial_index+1, valstep=1)

    # スキャン条件表示用テキストを初期化
    ax.text(0.4, 1.05, f'ℏω: {QE_sets[initial_index][0]} meV, h: {QE_sets[initial_index][1]}, k: {QE_sets[initial_index][2]}, l: {QE_sets[initial_index][3]}', 
                                   horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    
    def update(val):
        # スライダーの値に基づいて角度セットを選択
        index = int(val)-1
        
        ax.clear()
        
        # 格子点を描画
        ax.scatter(grid_points[:, 0], grid_points[:, 1], color='black', s=10, label='Lattice Points')
        
        if fixe==0: # ei fix
            Ei = bpe
            Ef = bpe - QE_sets[index][0]
        elif fixe==1: # ef fix
            Ei = bpe + QE_sets[index][0]
            Ef = bpe
        ki_cal=(Ei/2.072)**(1/2)
        kf_cal=(Ef/2.072)**(1/2)
        
        # hkl_calベクトルの計算 (sp1v_2d と sp2v_2d の線形和)
        hkl_cal = QE_sets[index][1] * astar + QE_sets[index][2] * bstar + QE_sets[index][3] * cstar
        
        # 線形結合を解く
        c1, c2 = solve_linear_combination(sp1v, sp2v, hkl_cal)
        
        hkl_x = c1 * sp1v_2d[0] + c2 * sp2v_2d[0]
        hkl_y = c1 * sp1v_2d[1] + c2 * sp2v_2d[1]
        # hkl_calベクトルを描画
        ax.quiver(0, 0, hkl_x, hkl_y, angles='xy', scale_units='xy', scale=1, color='red', label='τ (ki-kf)')
        
        # sp1v_2d, sp2v_2dベクトルを描画
        # 端に表示するようにする。
        ax.quiver(-sp1v_2d[0]*n_points1-sp2v_2d[0]*n_points2, -sp1v_2d[1]*n_points1-sp2v_2d[1]*n_points2, sp1v_2d[0], sp1v_2d[1], angles='xy', scale_units='xy', scale=1, color='magenta', label='axis 1')
        ax.quiver(-sp1v_2d[0]*n_points1-sp2v_2d[0]*n_points2, -sp1v_2d[1]*n_points1-sp2v_2d[1]*n_points2, sp2v_2d[0], sp2v_2d[1], angles='xy', scale_units='xy', scale=1, color='aqua', label='axis 2')
        
        # ki, kfベクトルを描画. ki//y_axis
        #C_sets.append([C1, C2, C3, C4])
        #inst_x = ki_cal * np.sin(np.radians(C_sets[1]-C_sets[3]))
        #inst_y = ki_cal * np.cos(np.radians(C_sets[1]-C_sets[3]))
        ki_vx = ki_cal * np.sin(np.radians(C_sets[index][1]-C_sets[index][3]))
        ki_vy = ki_cal * np.cos(np.radians(C_sets[index][1]-C_sets[index][3]))
        
        kf_vx = kf_cal * np.sin(np.radians(C_sets[index][1]-C_sets[index][3]+A_sets[index][1]))
        kf_vy = kf_cal * np.cos(np.radians(C_sets[index][1]-C_sets[index][3]+A_sets[index][1]))
        ax.quiver(hkl_x-ki_vx, hkl_y-ki_vy, ki_vx, ki_vy,  angles='xy', scale_units='xy', scale=1, color='blue', label='ki')
        ax.quiver(hkl_x-ki_vx, hkl_y-ki_vy, kf_vx, kf_vy, angles='xy', scale_units='xy', scale=1, color='green', label='kf')
        
        # 軸の設定
        #plt.xlim(-1.5 * max(norm_sp1v, norm_sp2v), 1.5 * max(norm_sp1v, norm_sp2v))
        #plt.ylim(-1.5 * max(norm_sp1v, norm_sp2v), 1.5 * max(norm_sp1v, norm_sp2v))
        
        ax.set_aspect('equal', adjustable='box')
        
        # ラベルとタイトル
        #ax.title('Reciprocal Space Plot')
        ax.set_xlabel(r'$k_x\ (\mathrm{\AA}^{-1})$')
        ax.set_ylabel(r'$k_y\ (\mathrm{\AA}^{-1})$')
        
        # スキャン条件表示用テキストを初期化
        ax.text(0.4, 1.05, f'ℏω: {QE_sets[index][0]} meV, h: {QE_sets[index][1]}, k: {QE_sets[index][2]}, l: {QE_sets[index][3]}', 
                                    horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

        
        # 凡例をグラフの外に配置
        ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0., fontsize=10)
        
        plt.draw()

    slider.on_changed(update)

    # 初期の描画を行う
    update(slider.val)

    def on_key(event):
        """
        左右矢印キーでスライダーを動かす
        """
        if event.key == 'right':
            slider.set_val(min(slider.val + 1, len(QE_sets)))
        elif event.key == 'left':
            slider.set_val(max(slider.val - 1, 1))

    # キーイベントを設定
    fig.canvas.mpl_connect('key_press_event', on_key)
    
    #plt.show()
