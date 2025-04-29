import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import pi, sin, cos, tan, sqrt, arcsin, arccos
from numpy import arctan2  # atan2はarctan2としてインポートする必要があります
import configparser
import os
import sys
from matplotlib.widgets import Slider
from scipy.optimize import minimize_scalar
from scipy.optimize import minimize
import pandas as pd

from PIL import Image  # GIF 保存のために必要

def calcresolution_scan(A_sets,QE_sets,Ni_mir,bpe,fixe,Hfocus,num_ana,entry_values,initial_index=0,save_gif=False,gif_name="resolution.gif"):
    # save_gifがTrueだと保存、falseだと非保存

    # INIファイルから設定を読み込む
    config = configparser.ConfigParser()
    # .exe化した場合に対応する
    if getattr(sys, 'frozen', False):
        # .exeの場合、sys.argv[0]が実行ファイルのパスになる
        ini_path = os.path.join(os.path.dirname(sys.argv[0]), 'config.ini')
    else:
        # .pyの場合、__file__がスクリプトのパスになる
        ini_path = os.path.join(os.path.dirname(__file__), 'config.ini')

    config.read(ini_path)
    
    view_mode = config['settings']['system']
    
    # divergenceの読み出し
    div_1st_m = float(entry_values.get("div_1st_m"))
    div_1st_h = float(entry_values.get("div_1st_h"))
    div_1st_v = float(entry_values.get("div_1st_v"))
    div_2nd_h = float(entry_values.get("div_2nd_h"))
    div_2nd_v = float(entry_values.get("div_2nd_v"))
    div_3rd_h = float(entry_values.get("div_3rd_h"))
    div_3rd_v = float(entry_values.get("div_3rd_v"))
    div_4th_h = float(entry_values.get("div_4th_h"))
    div_4th_v = float(entry_values.get("div_4th_v"))
    
    mos_mono_h = float(entry_values.get("mos_mono_h"))
    mos_mono_v = float(entry_values.get("mos_mono_v"))
    mos_sam_h = float(entry_values.get("mos_sam_h"))
    mos_sam_v = float(entry_values.get("mos_sam_v"))
    mos_ana_h = float(entry_values.get("mos_ana_h"))
    mos_ana_v = float(entry_values.get("mos_ana_v"))
    
    sample_to_analyzer = float(config['settings']['sample_to_analyzer'])
    analyzer_width = float(config['settings']['analyzer_width'])
    
    # プロット設定
    # グラフの描画
    fig, ax = plt.subplots(figsize=(8, 5))
    plt.subplots_adjust(left=0.15, bottom=0.25)  # 左の余白を削る

    # スライダー設定
    ax_slider = plt.axes([0.25, 0.10, 0.65, 0.03], facecolor='lightgoldenrodyellow')
    slider = Slider(ax_slider, 'scan number', 1, len(A_sets) , valinit=initial_index+1, valstep=1)
    
    # フレーム保存用リスト
    frames = []
    
    # グラフのタイトル（楕円の説明）
    ax.set_title("red circle : $Q_{\\parallel}$, blue circle : $Q_{\\perp}$", fontsize=12)

    # 追加情報を ax.text で追加
    ax.text(
        0.5, 1.1,  # グラフの外に配置 (x=0.4, y=1.05)
        f'ℏω: {QE_sets[initial_index][0]} meV, h: {QE_sets[initial_index][1]}, k: {QE_sets[initial_index][2]}, l: {QE_sets[initial_index][3]}',
        horizontalalignment='center',
        verticalalignment='center',
        transform=ax.transAxes,
        fontsize=10
    )
    
    def update(val):
        index = int(val)-1
        A1, A2, A3 = A_sets[index]
        hw = QE_sets[index][0]
        
        
        # プロットを再描画
        ax.clear()
        
        # system設定に基づいて角度の係数を変更
        EM = 1
        EA = 1 #anti-W配置の場合-1 

        if fixe==0: # ei fix
            Ei = bpe
            Ef = bpe - hw
        elif fixe==1: # ef fix
            Ei = bpe + hw
            Ef = bpe

        ki=(Ei/2.072)**(1/2)
        kf=(Ef/2.072)**(1/2)

        Q = np.sqrt(ki**2 + kf**2 - 2 * ki * kf * np.cos(np.radians(A2))) 

        # C1とA1の計算
        C1 = A1/2

        # C3とA3の計算
        C3 = A3/2

        thetaM = C1 * EM
        thetaS = A2 / 2
        thetaA = C3 * EA
        
        phi = np.degrees(np.arctan2(-kf * np.sin(np.radians(2 * thetaS)), ki - kf * np.cos(np.radians(2 * thetaS))))

        # Define constants for the resolution matrices
        # ここでαi とβi は、コリメータの水平方向と鉛直方向における発散角を表している。η とη′をモノクロメータとアナライザの水平方向と鉛直方向のモザイクのFWHM

        #theta0 = 0.1 #(A^-1)
        #lamda = (81.81 / Ei)**(1/2)
        # 0.4246609 = 1/(2*sqrt(2*log(2)))
        if Ni_mir == 0:
            alpha1 = div_1st_h / 60 / 180 * pi * 0.4246609
            beta1 = div_1st_v / 60 / 180 * pi * 0.4246609 
        elif Ni_mir == 1:
            alpha1 = div_1st_m*2*arcsin(0.0219*(81.81/Ei)**(1/2)/(4*pi))/pi*180*60*0.4246609 # NiのQcは0.0219
            beta1 = alpha1
            
        alpha2 = div_2nd_h / 60 / 180 * pi * 0.4246609
        # focusingの場合式が異なる。
        if Hfocus==0:
            alpha3 = div_3rd_h / 60 / 180 * pi * 0.4246609
        elif Hfocus==1:
            L=sample_to_analyzer
            W=analyzer_width*num_ana*np.sin(np.radians(A3))
            af=2 * np.degrees(np.arctan((W / 2) / L))
            #alpha3 = div_3rd_h / 60 / 180 * pi * 0.4246609 * (8*np.log(2)/12)**(1/2)
            alpha3 = af / 180 * pi * 0.4246609 * (8*np.log(2)/12)**(1/2)
        
        alpha4 = div_4th_h / 60 / 180 * pi * 0.4246609
        beta2 = div_2nd_v / 60 / 180 * pi * 0.4246609
        beta3 = div_3rd_v / 60 / 180 * pi * 0.4246609
        beta4 = div_4th_v / 60 / 180 * pi * 0.4246609
        
        etaM = mos_mono_h / 60 / 180 * pi * 0.4246609
        etaA = mos_ana_h / 60 / 180 * pi * 0.4246609
        etaS = mos_sam_h / 60 / 180 * pi * 0.4246609
        etaMp = mos_mono_v / 60 / 180 * pi * 0.4246609
        etaAp = mos_ana_v / 60 / 180 * pi * 0.4246609
        etaSp = mos_sam_v / 60 / 180 * pi * 0.4246609

        # Gについてreslibと同じ値になることを確認
        G = 8 * np.log(2) / (8 * np.log(2)) * np.diag([1 / (alpha1 ** 2), 1 / (alpha2 ** 2), 1 / (beta1 ** 2), 1 / (beta2 ** 2), 1 / (alpha3 ** 2), 1 / (alpha4 ** 2), 1 / (beta3 ** 2), 1 / (beta4 ** 2)])
        # Fについてreslibと同じ値になることを確認
        F = 8 * np.log(2) / (8 * np.log(2))  * np.diag([1 / (etaM ** 2), 1 / (etaMp ** 2), 1 / (etaA ** 2), 1 / (etaAp ** 2)])
        
        # Define matrices A, B, and C
        A = np.zeros((6, 8))# reslibと一致することを確認
        C = np.zeros((4, 8))# reslibと一致することを確認
        B = np.zeros((4, 6))# reslibと一致することを確認
        
        A[0, 0] = ki / (2 * np.tan(np.radians(thetaM)))
        A[0, 1] = -A[0, 0]
        A[3, 4] = kf / (2 * np.tan(np.radians(thetaA)))
        A[3, 5] = -A[3, 4]
        A[1, 1] = ki
        A[2, 3] = ki
        A[4, 4] = kf
        A[5, 6] = kf
        
        # 2.072142=h^2/m
        B[0, 0] = np.cos(np.radians(phi))
        B[0, 1] = np.sin(np.radians(phi))
        B[0, 3] = -np.cos(np.radians(phi - 2 * thetaS))
        B[0, 4] = -np.sin(np.radians(phi - 2 * thetaS))
        B[1, 0] = -B[0, 1]
        B[1, 1] = B[0, 0]
        B[1, 3] = -B[0, 4]
        B[1, 4] = B[0, 3]
        B[2, 2] = 1
        B[2, 5] = -1
        B[3, 0] = 2 * 2.072142 * ki
        B[3, 3] = -2 * 2.072142 * kf

        C[0, 0] = 1 / 2
        C[0, 1] = 1 / 2
        C[2, 4] = 1 / 2
        C[2, 5] = 1 / 2
        C[1, 2] = 1 / (2 * sin(np.radians(thetaM)))
        C[1, 3] = -C[1, 2]
        C[3, 6] = 1 / (2 * sin(np.radians(thetaA)))
        C[3, 7] = -C[3, 6]
        
        # 計算
        term = np.linalg.inv(G + C.T @ F @ C)  # G + C' * F * C の逆行列
        HF = A @ term @ A.T  # A * (G + C' * F * C)^(-1) * A'
        # HFまでreslibと一致
        if Hfocus == 1:
            P = np.linalg.inv(HF)
            P[4, 4] = (1 / (kf * alpha3)) ** 2
            P[3, 4] = 0
            P[3, 3] = (np.tan(np.radians(thetaA)) / (etaA * kf)) ** 2
            P[4, 3] = 0
            Pinv = np.linalg.inv(P)
            Minv = B @ Pinv @ B.T
        if Hfocus == 0:
            Minv = B @ HF @ B.T #これもreslibと一致
        M = np.linalg.inv(Minv)
        
        # RM 行列の設定
        #RM = np.zeros((4, 4))  # 4x4 のゼロ行列で初期化
        """
        RM = [[M[0, 0],M[0, 1],M[0, 3],M[0, 2]],
            [M[1, 0],M[1, 1],M[1, 3],M[1, 2]],
            [M[3, 0],M[3, 1],M[3, 3],M[3, 2]],
            [M[2, 0],M[2, 1],M[2, 3],M[2, 2]]]
        """
        # 軸 2↔3 をスワップするインデックス
        swap = [0, 1, 3, 2]
        RM = M[np.ix_(swap, swap)]
        
        # サンプルモザイクを入れた場合の計算
        Minv = np.linalg.inv(RM)
        Minv[1, 1] += Q**2 * etaS**2# / (8 * np.log(2))
        Minv[3, 3] += Q**2 * etaSp**2# / (8 * np.log(2))
        RM = np.linalg.inv(Minv)
        
        # RMは(q//,q⊥,hw,qz)における空間分布
        """
        # 保存するリスト
        data_list = []

        # ヘッダーと行列を順に追加する関数
        def add_matrix_to_list(name, matrix):
            data_list.append([f"=== {name} ==="])
            for row in matrix:
                data_list.append(list(row))
            data_list.append([])  # 空行を挿入

        # 追加していく
        add_matrix_to_list('A', A)
        add_matrix_to_list('B', B)
        add_matrix_to_list('C', C)
        #add_matrix_to_list('D', D)
        add_matrix_to_list('G', G)
        add_matrix_to_list('F', F)
        add_matrix_to_list('HF', HF)
        add_matrix_to_list('M', M)
        add_matrix_to_list('Minv', Minv)
        add_matrix_to_list('RM', RM)

        # DataFrameに変換（短い行にはNaNが入る）
        df = pd.DataFrame(data_list)

        # CSV出力
        df.to_csv('matrices_output.csv', index=False, header=False)
        """
        # プロット範囲
        #Xrange_lim = 0.1
        #Zrange_lim = 0.5
        
        # Qx=Q//,Qy=Q⊥の定義
        
        # 投影図の楕円の係数を計算する関数
        # fun3=@(x,y,z) RM(1,1).*x.^2+RM(2,2).*y.^2+RM(3,3).*z.^2+2*RM(1,2).*x.*y+2*RM(1,3).*x.*z+2*RM(2,3).*y.*z-2*log(2);
        def ellipse_coefficients(RM, log2, plane="xz"):
            if plane == "xz":
                A = RM[0, 0]
                C = RM[2, 2]
                B = 2 * RM[0, 2]
                D = 0  # xの線形項
                E = 0  # zの線形項
            elif plane == "yz":
                A = RM[1, 1]
                C = RM[2, 2]
                B = 2 * RM[1, 2]
                D = 0  # yの線形項
                E = 0  # zの線形項


            F = -2 * log2

            # 消去された軸（影響を除外）を考慮
            if plane == "xz":
                correction = RM[1, 1]  # y軸の効果を除外
                if correction != 0:
                    A -= (RM[0, 1]**2) / correction
                    C -= (RM[1, 2]**2) / correction
                    B -= 2 * (RM[0, 1] * RM[1, 2]) / correction
            elif plane == "yz":
                correction = RM[0, 0]  # x軸の効果を除外
                if correction != 0:
                    A -= (RM[0, 1]**2) / correction
                    C -= (RM[0, 2]**2) / correction
                    B -= 2 * (RM[0, 1] * RM[0, 2]) / correction
            return A, B, C, D, E, F

        # 楕円をプロットする関数
        def plot_ellipse(A, B, C, D, E, F, Xrange_lim, Zrange_lim, label, color,shift_x=0,shift_y=0):
            x = np.linspace(-Xrange_lim, Xrange_lim, 500)
            z = np.linspace(-Zrange_lim, Zrange_lim, 500)
            X, Z = np.meshgrid(x, z)

            # 楕円の式
            ellipse = A * X**2 + B * X * Z + C * Z**2 + D * X + E * Z + F
            
            # y方向にhwだけずらす
            X_shifted = X + shift_x
            # y方向にhwだけずらす
            Z_shifted = Z + shift_y

            # 等高線をプロット（楕円の曲線部分）
            #plt.contour(X_shifted, Z_shifted, ellipse, levels=[0], colors=color, label=label)
            ax.contour(X_shifted, Z_shifted, ellipse, levels=[0], colors=color, label=label)

        log2 = np.log(2)

        # xz平面の楕円の係数
        A_xz, B_xz, C_xz, D_xz, E_xz, F_xz = ellipse_coefficients(RM, log2, plane="xz")

        # yz平面の楕円の係数
        A_yz, B_yz, C_yz, D_yz, E_yz, F_yz = ellipse_coefficients(RM, log2, plane="yz")

        # 楕円球の係数行列 RM と楕円球の方程式
        def fun3(x, y, z, RM):
            return (
                RM[0, 0] * x**2
                + RM[1, 1] * y**2
                + RM[2, 2] * z**2
                + 2 * RM[0, 1] * x * y
                + 2 * RM[0, 2] * x * z
                + 2 * RM[1, 2] * y * z
                - 2 * np.log(2)
            )
        
        # 制約条件（楕円球の式 = 0 を満たす）
        def constraint(params, RM):
            x, y, z = params
            return fun3(x, y, z, RM)
        
        # 最大値を探索する関数
        def find_max_along_axis(RM, axis="x"):
            # 初期値
            initial_guess = [0, 0, 0]  # 楕円球の中心に近い点から探索を開始
            axis_map = {"x": 0, "y": 1, "z": 2}
            idx = axis_map[axis]

            # 目的関数（探索する軸を最大化）
            def objective(params):
                return -params[idx]  # 最大化したいのでマイナスを付ける

            # 制約条件を定義
            constraints = {"type": "eq", "fun": constraint, "args": (RM,)}

            # 最適化
            result = minimize(
                objective,
                initial_guess,
                method="SLSQP",
                constraints=constraints,
                options={"disp": False},
            )
            return result.x[idx], result.x  # 最大値とそのときの座標
        
        # 各軸の最大値を計算
        max_x, coords_x = find_max_along_axis(RM, axis="x")# Q//
        max_y, coords_y = find_max_along_axis(RM, axis="y")# Q⊥
        max_z, coords_z = find_max_along_axis(RM, axis="z")# E
        
        if max_y>max_x:
            Xrange_lim=max_y*1.5
        elif max_y<max_x:
            Xrange_lim=max_x*1.5
            
        Zrange_lim=max_z*1.5
        
        # xz平面とyz平面の楕円を描画
        plot_ellipse(A_xz, B_xz, C_xz, D_xz, E_xz, F_xz, Xrange_lim, Zrange_lim, label = "", color="red",shift_x=0, shift_y=0)
        plot_ellipse(A_yz, B_yz, C_yz, D_yz, E_yz, F_yz, Xrange_lim, Zrange_lim, label = "", color="blue",shift_x=0, shift_y=0)
        
        # 各軸の最大値を2倍した値
        resolution_Q_parallel = 2 * max_x
        resolution_Q_perpendicular = 2 * max_y
        resolution_energy = 2 * max_z
        
        # 軸やラベルの設定
        ax.axhline(0, color="black", linestyle="--", linewidth=0.5)
        ax.axvline(0, color="black", linestyle="--", linewidth=0.5)
        ax.set_xlabel("$Q$ ($\AA^{-1}$)")
        ax.set_ylabel("ℏω (meV)")
        
        # グラフのタイトル（楕円の説明）
        ax.set_title("red circle : $Q_{\\parallel}$, blue circle : $Q_{\\perp}$", fontsize=12)
        #ax.set_title("red circle : $Q_{x}$, blue circle : $Q_{y}$", fontsize=12)

        # 追加情報を ax.text で追加
        ax.text(
            0.5, 1.1,  # グラフの外に配置 (x=0.4, y=1.05)
            f'ℏω: {QE_sets[index][0]} meV, h: {QE_sets[index][1]}, k: {QE_sets[index][2]}, l: {QE_sets[index][3]}, δ$Q_{{\\parallel}}$ = {resolution_Q_parallel:.4f}, δ$Q_{{\\perp}}$ = {resolution_Q_perpendicular:.4f}, δE = {resolution_energy:.4f}',
            #f'ℏω: {QE_sets[index][0]} meV, h: {QE_sets[index][1]}, k: {QE_sets[index][2]}, l: {QE_sets[index][3]}, δ$Q_{{x}}$ = {resolution_Q_parallel:.4f}, δ$Q_{{y}}$ = {resolution_Q_perpendicular:.4f}, δE = {resolution_energy:.4f}',
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes,
            fontsize=10
        )

        # 軸の表示範囲を明示的に設定
        ax.set_xlim([-Xrange_lim, Xrange_lim])
        ax.set_ylim([-Zrange_lim, Zrange_lim])
        ax.grid(True)
        
        # フレーム保存（GIF用）
        if save_gif:
            fig.canvas.draw()
            image = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
            image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
            frames.append(Image.fromarray(image))

        # プロットの表示
        plt.draw()
        
    slider.on_changed(update)

    # 初期の描画を行う
    update(slider.val)
    
    def on_key(event):
        """
        左右矢印キーでスライダーを動かす
        """
        if event.key == 'right':
            slider.set_val(min(slider.val + 1, len(A_sets)))
        elif event.key == 'left':
            slider.set_val(max(slider.val - 1, 1))

    # キーイベントを設定
    fig.canvas.mpl_connect('key_press_event', on_key)
    
    if save_gif:
        for val in range(1, len(A_sets) + 1):
            slider.set_val(val)
        frames[0].save(gif_name, save_all=True, append_images=frames[1:], duration=100, loop=0)
        print(f"GIF 保存完了: {gif_name}")