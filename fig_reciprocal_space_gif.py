import configparser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d import Axes3D
import os
import sys
from PIL import Image  # GIF 保存のために必要

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

def plot_reciprocal_space_with_gif(bpe, bpc2, cphw, cp, fixe, sv1, sv2, RLtable,A_sets, C_sets,QE_sets, initial_index=0,save_gif=True,gif_name="reciprocalspace.gif"):
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
    hkl_cal = QE_sets[initial_index][1] * astar + QE_sets[initial_index][2] * bstar + QE_sets[initial_index][3] * cstar
    
    # 線形結合を解く
    c1, c2 = solve_linear_combination(sp1v, sp2v, hkl_cal)
    
    # プロットの準備
    fig, ax = plt.subplots(figsize=(8, 5))  # 図のサイズを縮小
    plt.subplots_adjust(left=0.01, bottom=0.25)
    
    # スライダー設定
    ax_slider = plt.axes([0.25, 0.10, 0.65, 0.03], facecolor='lightgoldenrodyellow')
    slider = Slider(ax_slider, 'scan number', 1, len(QE_sets) , valinit=initial_index+1, valstep=1)
    
    # フレーム保存用リスト
    frames = []

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
        
        # フレーム保存（GIF用）
        if save_gif:
            fig.canvas.draw()
            image = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
            image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
            frames.append(Image.fromarray(image))

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
    
    if save_gif:
        for val in range(1, len(A_sets) + 1):
            slider.set_val(val)
        frames[0].save(gif_name, save_all=True, append_images=frames[1:], duration=100, loop=0)
        print(f"GIF 保存完了: {gif_name}")
    
    #plt.show()
