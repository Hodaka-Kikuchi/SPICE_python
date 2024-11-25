import math
import numpy as np
from scipy.optimize import minimize
import configparser
import os
import sys

def angle_calc2(astar,bstar,cstar,UB,bpe,bpc2,bpmu,bpnu,bp,fixe,hw_ini,hw_fin,hw_inc,h_cal,k_cal,l_cal):
    # bragg peakの位置からoffsetを算出
    hkl_bp=bp[0]*astar+bp[1]*bstar+bp[2]*cstar
    #計算されたrlu
    N_hkl_bp=np.linalg.norm(hkl_bp)
    #A2を計算
    kf=(bpe/2.072)**(1/2)
    ki=(bpe/2.072)**(1/2)
    phi_bp = np.degrees(np.arccos((ki**2 + kf**2 - N_hkl_bp**2) / (2 * ki * kf)))
    # thetaの計算.Qとkiのなす角度
    theta_bp = np.degrees(np.arctan((ki - kf * np.cos(np.radians(phi_bp))) / (kf * np.sin(np.radians(phi_bp)))))
    Theta_bp_mat =  np.array([[np.cos(np.radians(theta_bp)), -np.sin(np.radians(theta_bp)), 0],[np.sin(np.radians(theta_bp)),  np.cos(np.radians(theta_bp)), 0],[0, 0, 1]])
    # QLの計算
    QL_bp = np.array([0, kf, 0]) - np.array([-kf * np.sin(np.radians(phi_bp)), kf * np.cos(np.radians(phi_bp)), 0])
    # Qthetaの計算
    Qtheta_bp = np.linalg.inv(Theta_bp_mat)@QL_bp
    Qtheta_bp = Qtheta_bp.reshape(-1, 1)  # 列ベクトルに変
    Qtheta_bp[np.abs(Qtheta_bp) <= 1e-6] = 0 #超重要,他のものにも適応
    # UBの更新
    # tiltの情報もここに入れてしまう。SPICEと同じ定義。
    N_bp=np.array([[1,0,0],[0,np.cos(np.radians(bpnu)),-np.sin(np.radians(bpnu))],[0,np.sin(np.radians(bpnu)),np.cos(np.radians(bpnu))]])
    M_bp=np.array([[np.cos(np.radians(bpmu)),0,np.sin(np.radians(bpmu))],[0,1,0],[-np.sin(np.radians(bpmu)),0,np.cos(np.radians(bpmu))]])
    UBt = 2 * np.pi * UB
    # Qv_pbの計算
    Qv_bp=UBt@(np.array([bp[0], bp[1], bp[2]]))
    Qv_bp = Qv_bp.reshape(-1, 1)  # 列ベクトルに変
    Qv_bp[np.abs(Qv_bp) <= 1e-6] = 0 #超重要,他のものにも適応
    
    # fitting process
    def Omera_toration(omega):
        return np.array([[np.cos(np.radians(omega)),-np.sin(np.radians(omega)),0],[np.sin(np.radians(omega)),np.cos(np.radians(omega)),0],[0,0,1]])
    def N_rotation(nu):
        return np.array([[1,0,0],[0,np.cos(np.radians(nu)),-np.sin(np.radians(nu))],[0,np.sin(np.radians(nu)),np.cos(np.radians(nu))]])
    def M_rotation(mu):
        return np.array([[np.cos(np.radians(mu)),0,np.sin(np.radians(mu))],[0,1,0],[-np.sin(np.radians(mu)),0,np.cos(np.radians(mu))]])
    
    # 最適化対象の関数
    def objective_bp(angles, Qv_bp, Qtheta_bp):
        omega, mu, nu = angles
        rotation_matrix0 = Omera_toration(omega) @ N_rotation(mu) @ M_rotation(nu)
        transformed_vector0 = rotation_matrix0 @ Qv_bp
        return np.linalg.norm(transformed_vector0 - Qtheta_bp)

    # 最適化関数（omegaの各初期値で最適化を行い、最も低い誤差を選択）
    def optimize_with_fixed_omega_bp(Qv_bp, Qtheta_bp, omega_values):
        best_result_bp = None
        best_fun_bp = float('inf')  # 最小値を探す

        for omega in omega_values:
            initial_guess = [omega, bpmu, bpnu]  # omegaを固定し、mu, nuは初期値0で開始
            result_bp = minimize(objective_bp, initial_guess, method='L-BFGS-B', bounds=[(-180, 180), (-90, 90), (-90, 90)], args=(Qv_bp, Qtheta_bp))
            
            # 最小の最適化結果を選ぶ
            if result_bp.fun < best_fun_bp:
                best_result_bp = result_bp
                best_fun_bp = result_bp.fun

        return best_result_bp

    # omegaの候補値
    omega_values = [-180, -135, -90, -45, 0, 45, 90, 135]
    
    # 最適化結果を格納している変数
    best_result_bp = optimize_with_fixed_omega_bp(Qv_bp, Qtheta_bp, omega_values)
    
    # 結果を表示
    omega0 = best_result_bp.x[0]
    
    # 結果を表示(-180~180に規格化。)
    #omega0, mu0, nu0 = [(angle + 180) % 360 - 180 for angle in result0.x]
    
    # s_iniの計算
    s_ini = omega0 + theta_bp
    # offsetの計算
    offset = bpc2 - s_ini
    
    ####################################################################################################
    # 変数はhw
    hw_tab = np.arange(hw_ini, hw_fin, hw_inc)
    if hw_tab[-1] != hw_fin:  # 最後の点が hw_fin でない場合
        hw_tab = np.append(hw_tab, hw_fin)

    
    # hklの位置から各angleを算出
    hkl_cal=h_cal*astar+k_cal*bstar+l_cal*cstar
    Nhkl_cal=np.linalg.norm(hkl_cal)
    
    results = []  # 結果を保存するためのリスト
    
    for i in range(len(hw_tab)):
        if fixe==0: # ei fix
            Ei = bpe
            Ef = bpe - hw_tab[i]
        elif fixe==1: # ef fix
            Ei = bpe + hw_tab[i]
            Ef = bpe
        #phi_cal = np.degrees(np.arccos((ki**2 + kf_cal**2 - Nhkl_cal**2) / (2 * ki_cal * kf_cal)))
        ki_cal=(Ei/2.072)**(1/2)
        kf_cal=(Ef/2.072)**(1/2)
        
        # phi_calの計算
        phi_cal = np.degrees(np.arccos((ki_cal**2 + kf_cal**2 - Nhkl_cal**2) / (2 * ki_cal * kf_cal)))
        
        # phi_calが正常に計算された場合のみ、後続の計算を実行
        theta_cal = np.degrees(np.arctan((ki_cal - kf_cal * np.cos(np.radians(phi_cal))) / (kf_cal * np.sin(np.radians(phi_cal)))))
        Theta_cal_mat =  np.array([[np.cos(np.radians(theta_cal)), -np.sin(np.radians(theta_cal)), 0],[np.sin(np.radians(theta_cal)),  np.cos(np.radians(theta_cal)), 0],[0, 0, 1]])
        QL_cal = np.array([0, ki_cal, 0]) - np.array([-kf_cal * np.sin(np.radians(phi_cal)), kf_cal * np.cos(np.radians(phi_cal)), 0])
        Qtheta_cal = np.linalg.inv(Theta_cal_mat)@(QL_cal)
        Qtheta_cal=Qtheta_cal.reshape(-1, 1)
        Qtheta_cal[np.abs(Qtheta_cal) <= 1e-6] = 0 #超重要,他のものにも適応
        Qv_cal = UBt@(np.array([h_cal, k_cal, l_cal]))
        Qv_cal=Qv_cal.reshape(-1, 1)
        Qv_cal[np.abs(Qv_cal) <= 1e-6] = 0 #超重要,他のものにも適応
    
        # 最適化対象の関数
        def objective(angles, Qv_cal, Qtheta_cal):
            omega, mu, nu = angles
            rotation_matrix = Omera_toration(omega) @ N_rotation(mu) @ M_rotation(nu)
            transformed_vector = rotation_matrix @ Qv_cal
            return np.linalg.norm(transformed_vector - Qtheta_cal)

        # 最適化関数（omegaの各初期値で最適化を行い、最も低い誤差を選択）
        def optimize_with_fixed_omega(Qv_cal, Qtheta_cal, omega_values):
            best_result = None
            best_fun = float('inf')  # 最小値を探す

            for omega in omega_values:
                initial_guess = [omega, 0, 0]  # omegaを固定し、mu, nuは初期値0で開始
                # argsにQv_calとQtheta_calを渡す
                result = minimize(
                    objective,
                    initial_guess,
                    method='L-BFGS-B',
                    bounds=[(-180, 180), (-90, 90), (-90, 90)],
                    args=(Qv_cal, Qtheta_cal)  # argsにQv_calとQtheta_calを渡す
                )

                # 最小の最適化結果を選ぶ
                if result.fun < best_fun:
                    best_result = result
                    best_fun = result.fun

            return best_result
        #最適化結果を格納している変数
        best_result = optimize_with_fixed_omega(Qv_cal, Qtheta_cal, omega_values)
        
        omega = best_result.x[0]
        mu = best_result.x[1]
        nu = best_result.x[2]
        
        s_cal=omega+theta_cal
        omega_inst=s_cal+offset
        if omega_inst>180:
            omega_inst = omega_inst-360
        elif omega_inst<-180:
            omega_inst = omega_inst+360
        
        # アナライザとモノクロメータの面間隔
        # INIファイルの設定読み込み
        config = configparser.ConfigParser()
        if getattr(sys, 'frozen', False):
            ini_path = os.path.join(os.path.dirname(sys.argv[0]), 'config.ini')
        else:
            ini_path = os.path.join(os.path.dirname(__file__), 'config.ini')
        config.read(ini_path)
        d_mono = float(config['instrument']['d_mono'])
        d_ana = float(config['instrument']['d_ana'])
        #d = 3.355  # PGの場合

        # C1とA1の計算
        C1 = np.degrees(np.arcsin((2 * np.pi / d_mono) / (2 * np.sqrt(Ei / 2.072))))
        A1 = 2 * C1
        
        # C3とA3の計算
        C3 = np.degrees(np.arcsin((2 * np.pi / d_ana) / (2 * np.sqrt(Ef / 2.072))))
        A3 = 2 * C3
        
        # 結果を辞書としてまとめる
        result = {
            'hw' : round(hw_tab[i],4),
            'h' : round(h_cal,4),
            'k' : round(k_cal,4),
            'l' : round(l_cal,4),
            'C1' : round(C1,4),
            'A1' : round(A1,4),
            'C2' : round(omega_inst,4),
            'A2' : round(phi_cal,4),
            'C3' : round(C3,4),
            'A3' : round(A3,4),
            'mu' : round(mu+bpmu,4),
            'nu' : round(nu+bpnu,4),
        }
        
        # リストに追加
        results.append(result)
    return results