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
from scipy.linalg import block_diag

def calcresolution_save(apr_value,sense,astar,bstar,cstar,sv1,sv2,A_sets,QE_sets,Ni_mir,bpe,fixe,focus_cond,entry_values):
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
    analyzer_width = float(config['settings']['ana_width'])

    # focusing conditionの読み出し
    MHF = focus_cond["MHF"]
    MVF = focus_cond["MVF"]
    AHF = focus_cond["AHF"]
    AVF = focus_cond["AVF"]

    num_mono_h = focus_cond["num_mono_h"]
    num_mono_v = focus_cond["num_mono_v"]
    num_ana_h  = focus_cond["num_ana_h"]
    num_ana_v  = focus_cond["num_ana_v"]

    reso_mat = np.zeros((4,4,len(A_sets)))
    col_cond = np.zeros((14,len(A_sets)))
    scan_cond = np.zeros((9,len(A_sets)))
    
    # matrix.shape[0]. shape はタプル (行数, 列数). shape[0] が行数、shape[1] が列数
    for index in range(len(A_sets)):
        A1, A2, A3 = A_sets[index]
        hw = QE_sets[index][0]
        
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
        #lamda = (81.81 / Ei)**(1/2) #(A^-1)
        # 0.4246609 = 1/(2*sqrt(2*log(2)))
        #alpha1 = div_1st_m * theta0 * lamda * ((2*np.log(2))**(1/2)) / (3*(1/2)) / 180 * pi * 0.4246609
        if Ni_mir == 0:
            alpha1 = div_1st_h / 60 / 180 * pi * 0.4246609
            beta1 = div_1st_v / 60 / 180 * pi * 0.4246609 
        elif Ni_mir == 1:
            NA = 6.022*10**(23) # mol^(-1)
            ro = 8.908 # g/cm^2
            M = 58.69 # g/mol
            bc = 1.03*10**(-12) # cm
            lamda=sqrt(81.81/Ei)*10**(-8) # cm
            Qc = 0.0219 # = sqrt(16*pi*po)
            div_1st_h = div_1st_m*2*np.degrees(arcsin(lamda*sqrt(NA*ro/M*bc/pi)))*60
            div_1st_v = div_1st_h
            alpha1 = div_1st_m*2*np.degrees(arcsin(lamda*sqrt(NA*ro/M*bc/pi)))*60  / 60 / 180 * pi * 0.4246609
            beta1 = alpha1
        alpha2 = div_2nd_h / 60 / 180 * pi * 0.4246609
        # focusingの場合式が異なる。
        # focusingの場合式が異なる。
        if AHF==0:
            alpha3 = div_3rd_h / 60 / 180 * pi * 0.4246609
        elif AHF==1:
            '''
            L=sample_to_analyzer
            W=analyzer_width*num_ana_h*np.sin(np.radians(A3))
            af=2 * np.degrees(np.arctan((W / 2) / L))
            #alpha3 = div_3rd_h / 60 / 180 * pi * 0.4246609 * (8*np.log(2)/12)**(1/2)
            alpha3 = af / 180 * pi * 0.4246609 * (8*np.log(2)/12)**(1/2)
            '''
            alpha3 = div_3rd_h / 60 / 180 * pi * 0.4246609 * (8*np.log(2)/12)**(1/2)
        
        col_cond[:,index]=[div_1st_h,div_2nd_h,div_3rd_h,div_4th_h,div_1st_v,div_2nd_v,div_3rd_v,div_4th_v,mos_mono_h,mos_ana_h,mos_sam_h,mos_mono_v,mos_ana_v,mos_sam_v]
        
        scan_cond[:,index]=[A1,A2,A3,Ei,Ef,hw,QE_sets[index][1], QE_sets[index][2], QE_sets[index][3]]
        
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
        
        # popovic近似へ分岐
        if apr_value == "P":
            beam_width = float(config['settings']['beam_width'])
            beam_height = float(config['settings']['beam_height'])
            # ==== Beam shape ====
            beamw = (beam_width)**2
            beamh = (beam_height)**2
            bshape = np.diag([beamw, beamh])

            # ==== Mono shape ====
            mono_width = float(config['settings']['mono_width'])
            mono_height = float(config['settings']['mono_height'])
            mono_depth = float(config['settings']['mono_depth'])
            monow = (num_mono_h*mono_width)**2
            monoh = (num_mono_v*mono_height)**2
            monod = (mono_depth)**2
            mshape = np.diag([monod, monow, monoh])

            '''
            # ==== Monitor shape ====
            # only flux normalization
            monitorw = 1
            monitorh = 1
            # if you want:
            #monitorw = monitor.width**2
            #monitorh = monitor.height**2
            monitorshape = np.diag([monitorw, monitorh])
            '''

            # ==== Sample shape ====
            sshape = np.eye(3)
            psi = thetaS -phi
            rot = np.array([[cos(np.radians(psi)),sin(np.radians(psi)),0],
                [-sin(np.radians(psi)),cos(np.radians(psi)),0],
                [0,0,1]])
            sshape = rot@sshape@rot.T

            # ==== Analyzer shape ====
            ana_width = float(config['settings']['ana_width'])
            ana_height = float(config['settings']['ana_height'])
            ana_depth = float(config['settings']['ana_depth'])
            anaw = (num_ana_h*ana_width)**2
            anah = (num_ana_v*ana_height)**2
            anad = (ana_depth)**2
            ashape = np.diag([anad, anaw, anah])

            # ==== Detector shape ====
            det_width = float(config['settings']['det_width'])
            det_height = float(config['settings']['det_height'])
            detectorw = (det_width)**2
            detectorh = (det_height)**2
            dshape = np.diag([detectorw, detectorh])

            # ==== S matrix ====
            Sinv = block_diag(bshape, mshape, sshape, ashape, dshape)  # S^-1
            S = np.linalg.inv(Sinv)

            # ==== Distances ====
            L0 = float(config['settings']['source_to_monochromator'])
            L1 = float(config['settings']['monochromator_to_sample'])
            #L1mon = 1 # only flux normalization
            L2 = float(config['settings']['sample_to_analyzer'])
            L3 = float(config['settings']['analyzer_to_detector'])

            def focusing_curvature(L_1, L_2, theta):
                # 有効焦点距離。ただし、単位をmmからmに直す必要がある。
                f = 1.0 / (1.0/L_1 + 1.0/L_2)

                # θ をラジアンに変換
                theta_rad = np.radians(theta)

                # 曲率 R = 2*f*|sin(theta)|
                R = 2.0 * f * np.abs(np.sin(theta_rad))

                return R

            # focusing:
            if MVF == 0:
                monorv = 1e10
            elif MVF ==1:
                monorv = focusing_curvature(L0,L1,thetaM)
            if MHF == 0:
                monorh = 1e10
            elif MHF ==1:
                monorh = focusing_curvature(L0,L1,thetaM)
            if AVF == 0:
                anarv = 1e10
            elif AVF ==1:
                anarv = focusing_curvature(L2,L3,thetaA)
            if AHF == 0:
                anarh = 1e10
            elif AHF ==1:
                anarh = focusing_curvature(L2,L3,thetaA)

            # ==== T matrix ====
            T = np.zeros((4, 13))

            T[0, 0] = -1/(2*L0)
            T[0, 2] = np.cos(np.radians(thetaM))*(1/L1 - 1/L0)/2
            T[0, 3] = np.sin(np.radians(thetaM))*(1/L0 + 1/L1 - 2/(monorh*np.sin(np.radians(thetaM))))/2
            T[0, 5] = np.sin(np.radians(thetaS))/(2*L1)
            T[0, 6] = np.cos(np.radians(thetaS))/(2*L1)

            T[1, 1] = -1/(2*L0*np.sin(np.radians(thetaM)))
            T[1, 4] = (1/L0 + 1/L1 - 2*np.sin(np.radians(thetaM))/monorv)/(2*np.sin(np.radians(thetaM)))
            T[1, 7] = -1/(2*L1*np.sin(np.radians(thetaM)))

            T[2, 5] = np.sin(np.radians(thetaS))/(2*L2)
            T[2, 6] = -np.cos(np.radians(thetaS))/(2*L2)
            T[2, 8] = np.cos(np.radians(thetaA))*(1/L3 - 1/L2)/2
            T[2, 9] = np.sin(np.radians(thetaA))*(1/L2 + 1/L3 - 2/(anarh*np.sin(np.radians(thetaA))))/2
            T[2, 11] = 1/(2*L3)

            T[3, 7] = -1/(2*L2*np.sin(np.radians(thetaA)))
            T[3, 10] = (1/L2 + 1/L3 - 2*np.sin(np.radians(thetaA))/anarv)/(2*np.sin(np.radians(thetaA)))
            T[3, 12] = -1/(2*L3*np.sin(np.radians(thetaA)))

            # ==== D matrix ====
            D = np.zeros((8, 13))

            D[0, 0] = -1/L0
            D[0, 2] = -np.cos(np.radians(thetaM))/L0
            D[0, 3] = np.sin(np.radians(thetaM))/L0

            D[2, 1] = D[0, 0]
            D[2, 4] = -D[0, 0]

            D[1, 2] = np.cos(np.radians(thetaM))/L1
            D[1, 3] = np.sin(np.radians(thetaM))/L1
            D[1, 5] = np.sin(np.radians(thetaS))/L1
            D[1, 6] = np.cos(np.radians(thetaS))/L1

            D[3, 4] = -1/L1
            D[3, 7] = -D[3, 4]

            D[4, 5] = np.sin(np.radians(thetaS))/L2
            D[4, 6] = -np.cos(np.radians(thetaS))/L2
            D[4, 8] = -np.cos(np.radians(thetaA))/L2
            D[4, 9] = np.sin(np.radians(thetaA))/L2

            D[6, 7] = -1/L2
            D[6, 10] = -D[6, 7]

            D[5, 8] = np.cos(np.radians(thetaA))/L3
            D[5, 9] = np.sin(np.radians(thetaA))/L3
            D[5, 11] = 1/L3

            D[7, 10] = -D[5, 11]
            D[7, 12] = D[5, 11]

        # 計算
        term = np.linalg.inv(G + C.T @ F @ C)  # G + C' * F * C の逆行列
        HF = A @ term @ A.T  # A * (G + C' * F * C)^(-1) * A'
        # HFまでreslibと一致
        if apr_value == "P":
            Minv = (B @ A @ np.linalg.inv(np.linalg.inv(D @ np.linalg.inv(S + T.T @ F @ T) @ D.T) + G) @ A.T @ B.T)
        elif apr_value == "CN":
            if AHF == 1:
                P = np.linalg.inv(HF)
                P[4, 4] = (1 / (kf * alpha3)) ** 2
                P[3, 4] = 0
                P[3, 3] = (np.tan(np.radians(thetaA)) / (etaA * kf)) ** 2
                P[4, 3] = 0
                Pinv = np.linalg.inv(P)
                Minv = B @ Pinv @ B.T
            elif AHF == 0:
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
        
        Qx = sv1[0]*astar+sv1[1]*bstar+sv1[2]*cstar
        Qy = sv2[0]*astar+sv2[1]*bstar+sv2[2]*cstar
        Qvect = QE_sets[index][1]*astar+QE_sets[index][2]*bstar+QE_sets[index][3]*cstar
        
        # Q方向の単位ベクトル
        uq = Qvect / np.linalg.norm(Qvect)

        # MATLABの scalar(u,v) に対応 → 普通の内積でOK (すでに直交座標系に変換済みなので)
        xq = np.dot(sv1, uq) / np.linalg.norm(sv1)
        yq = np.dot(sv2, uq) / np.linalg.norm(sv2)

        # 回転角度
        theta_rad = np.arctan2(yq, xq)

        # 回転行列（MATLABの tmat と同じ形）
        rot_mat = np.array([
            [ np.cos(theta_rad), -np.sin(theta_rad), 0, 0],
            [ np.sin(theta_rad),  np.cos(theta_rad), 0, 0],
            [0,                  0,                  1, 0],
            [0,                  0,                  0, 1]
        ])
        
        # 相似変換
        RM = rot_mat @ RM @ rot_mat.T
        
        if sense == '-+-':
            # 上下反転
            # rightではそのまま、leftでaxis2をaxis1に対してミラーさせる。
            S = np.diag([1.0, -1.0, 1.0, 1.0])   # y軸のみ反転
            RM_flipped = S @ RM @ S.T
            
            RM = RM_flipped
        elif sense == '+-+':
            pass
        
        # RMは(q//,q⊥,hw,qz)における空間分布
        # これを(qx(axis1),qy(axis2),hw,qz)に置ける空間分布に変換する。
        
        reso_mat[:,:,index] = RM
    return reso_mat,col_cond,scan_cond