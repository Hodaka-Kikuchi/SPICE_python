import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import pi, sin, cos, tan, sqrt, arcsin, arccos
from numpy import arctan2  # atan2はarctan2としてインポートする必要があります
import configparser
import os
import sys
from mpl_toolkits.mplot3d import Axes3D

def calcresolution(A_sets,QE_sets,bpe,fixe,hw,Hfocus,num_ana,entry_values):

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
    
    # 水平集光
    # Hfocus = 1
    # num_ana = 7
    
    # 設定を変数に代入
    """
    # instrumentセクションのview設定を読み込む
    div_1st_h = float(config['instrument']['div_1st_h'])
    div_1st_v = float(config['instrument']['div_1st_v'])
    div_2nd_h = float(config['instrument']['div_2nd_h'])
    div_2nd_v = float(config['instrument']['div_2nd_v'])
    div_3rd_h = float(config['instrument']['div_3rd_h'])
    div_3rd_v = float(config['instrument']['div_3rd_v'])
    div_4th_h = float(config['instrument']['div_4th_h'])
    div_4th_v = float(config['instrument']['div_4th_v'])

    mos_mono_h = float(config['instrument']['mos_mono_h'])
    mos_mono_v = float(config['instrument']['mos_mono_v'])
    mos_ana_h = float(config['instrument']['mos_ana_h'])
    mos_ana_v = float(config['instrument']['mos_ana_v'])

    mos_sam_h = float(config['sample']['mos_sam_h'])
    mos_sam_v = float(config['sample']['mos_sam_v'])
    
    d_mono = float(config['instrument']['d_mono'])
    d_ana = float(config['instrument']['d_ana'])
    """
    view_mode = config['settings']['system']
    
    # divergenceの読み出し
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
    
    A1, A2, A3 = A_sets[0]
    
    A2 = -A2 # 分光器の図を書くためにA2はマイナスにして値を渡しているため
    
    # system設定に基づいてy軸の設定を変更
    if view_mode == 'left':
        EM = -1
    elif view_mode == 'right':
        EM = 1

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

    thetaM = -C1 * EM
    thetaS = A2 / 2
    thetaA = -C3 * EM
    
    phi = np.degrees(np.arctan2(-kf * np.sin(np.radians(2 * thetaS)), ki - kf * np.cos(np.radians(2 * thetaS))))

    # Define constants for the resolution matrices
    # ここでαi とβi は、コリメータの水平方向と鉛直方向における発散角を表している。η とη′をモノクロメータとアナライザの水平方向と鉛直方向のモザイクのFWHM

    alpha1 = div_1st_h / 60 / 180 * pi
    alpha2 = div_2nd_h / 60 / 180 * pi
    # focusingの場合式が異なる。
    if Hfocus==0:
        alpha3 = div_3rd_h / 60 / 180 * pi
    elif Hfocus==1:
        L=sample_to_analyzer
        W=analyzer_width*num_ana
        af=np.degrees(2 * np.arctan((W / 2) / L))
        alpha3 = (8*np.log(2)/12)**(1/2)*af / 180 * pi
    
    alpha4 = div_4th_h / 60 / 180 * pi
    beta1 = div_1st_v / 60 / 180 * pi
    beta2 = div_2nd_v / 60 / 180 * pi
    beta3 = div_3rd_v / 60 / 180 * pi
    beta4 = div_4th_v / 60 / 180 * pi
    
    etaM = mos_mono_h / 60 / 180 * pi
    etaA = mos_ana_h / 60 / 180 * pi
    etaS = mos_sam_h / 60 / 180 * pi
    etaMp = mos_mono_v / 60 / 180 * pi
    etaAp = mos_ana_v / 60 / 180 * pi
    etaSp = mos_sam_v / 60 / 180 * pi

    G = 8 * np.log(2) * np.diag([1 / alpha1 ** 2, 1 / alpha2 ** 2, 1 / beta1 ** 2, 1 / beta2 ** 2, 1 / alpha3 ** 2, 1 / alpha4 ** 2, 1 / beta3 ** 2, 1 / beta4 ** 2])

    F = 8 * np.log(2) * np.diag([1 / etaM ** 2, 1 / etaMp ** 2, 1 / etaA ** 2, 1 / etaAp ** 2])

    # Define matrices A, B, and C
    A = np.zeros((6, 8))
    C = np.zeros((4, 8))
    B = np.zeros((4, 6))

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
    
    # 行列 C の転置（C' in MATLAB）
    C_T = C.T
    # 計算
    term = np.linalg.inv(G + C_T @ F @ C)  # G + C' * F * C の逆行列
    HF = A @ term @ A.T  # A * (G + C' * F * C)^(-1) * A'
    if Hfocus == 1:
        HF = np.linalg.inv(HF)
        HF[4, 4] = (1 / (kf * alpha3)) ** 2
        HF[4, 3] = 0
        HF[3, 4] = 0
        HF[3, 3] = (np.tan(np.radians(thetaA)) / (etaA * kf)) ** 2
        HF = np.linalg.inv(HF)
    Minv = B @ HF @ B.T
    M = np.linalg.inv(Minv)
    
    # RM 行列の設定
    RM = np.zeros((4, 4))  # 4x4 のゼロ行列で初期化
    RM[0, 0] = M[0, 0]
    RM[1, 0] = M[1, 0]
    RM[0, 1] = M[0, 1]
    RM[1, 1] = M[1, 1]
    RM[0, 2] = M[0, 3]
    RM[1, 2] = M[1, 3]
    RM[2, 0] = M[3, 0]
    RM[2, 2] = M[3, 3]
    RM[2, 1] = M[3, 1]
    RM[1, 2] = M[1, 3]
    RM[0, 3] = M[0, 2]
    RM[3, 0] = M[2, 0]
    RM[3, 3] = M[2, 2]
    RM[3, 1] = M[2, 1]
    RM[1, 3] = M[1, 2]

    # サンプルモザイクを入れた場合の計算
    # Minv の再計算
    Minv = np.linalg.inv(RM)
    
    # Minv の要素を修正
    Minv[1, 1] += Q**2 * etaS**2 / (8 * np.log(2))
    Minv[3, 3] += Q**2 * etaSp**2 / (8 * np.log(2))
    
    # RM の再計算
    RM = np.linalg.inv(Minv)
    # RMは(q//,q⊥,hw,qz)における空間分布
    # プロット範囲
    #Xrange_lim = 0.1
    #Zrange_lim = 0.5
    if Hfocus==0:
        Xrange_lim=Q*5/100
    elif Hfocus==1:
        Xrange_lim=Q*5*num_ana/100
    Zrange_lim=Ei*10/100
    
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
        plt.contour(X_shifted, Z_shifted, ellipse, levels=[0], colors=color, label=label)

    log2 = 2 * np.log(2)

    # xz平面の楕円の係数
    A_xz, B_xz, C_xz, D_xz, E_xz, F_xz = ellipse_coefficients(RM, log2, plane="xz")

    # yz平面の楕円の係数
    A_yz, B_yz, C_yz, D_yz, E_yz, F_yz = ellipse_coefficients(RM, log2, plane="yz")

    # グラフの描画
    fig, ax = plt.subplots(figsize=(8, 6))
    """
    # x=0 の場合（点線で表示）
    contour_x0 = plt.contour(Y, Z, F_x0, levels=[0], colors="blue", linestyles="--", label="x = 0")
    plt.plot([], [], color="blue", label='$Q_{y}$ ($\AA^{-1}$)')  # 凡例用

    # y=0 の場合（点線で表示）
    contour_y0 = plt.contour(X, Z, F_y0, levels=[0], colors="red", linestyles="--", label="y = 0")
    plt.plot([], [], color="red", label='$Q_{x}$ ($\AA^{-1}$)')  # 凡例用
    """
    # xz平面とyz平面の楕円を描画
    plot_ellipse(A_yz, B_yz, C_yz, D_yz, E_yz, F_yz, Xrange_lim, Zrange_lim, label = "", color="blue",shift_x=0, shift_y=0)
    plot_ellipse(A_xz, B_xz, C_xz, D_xz, E_xz, F_xz, Xrange_lim, Zrange_lim, label = "", color="red",shift_x=0, shift_y=0)

    # 軸やラベルの設定
    ax.axhline(0, color="black", linestyle="--", linewidth=0.5)
    ax.axvline(0, color="black", linestyle="--", linewidth=0.5)
    ax.set_xlabel("$Q$ ($\AA^{-1}$)")
    ax.set_ylabel("ℏω (meV)")
    
    # グラフのタイトル（楕円の説明）
    ax.set_title("red circle : $Q_{\\parallel}$, blue circle : $Q_{\\perp}$", fontsize=12)

    # 追加情報を ax.text で追加
    ax.text(
        0.5, 1.1,  # グラフの外に配置 (x=0.4, y=1.05)
        f'ℏω: {QE_sets[0][0]} meV, h: {QE_sets[0][1]}, k: {QE_sets[0][2]}, l: {QE_sets[0][3]}',
        horizontalalignment='center',
        verticalalignment='center',
        transform=ax.transAxes,
        fontsize=10
    )

    # 軸の表示範囲を明示的に設定
    ax.set_xlim([-Xrange_lim, Xrange_lim])
    ax.set_ylim([-Zrange_lim, Zrange_lim])
    ax.grid(True)

    # プロットの表示
    #plt.show()