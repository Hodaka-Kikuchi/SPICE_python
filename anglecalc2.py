import math
import numpy as np
from scipy.optimize import minimize

def angle_calc(astar,bstar,cstar,UB,bpe,bpc2,bpmu,bpnu,bp,cphw,cp,fixe,angle_ax1_ax2):
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
    # UBの更新
    # tiltの情報もここに入れてしまう。SPICEと同じ定義。
    N_bp=np.array([[1,0,0],[0,np.cos(np.radians(bpnu)),-np.sin(np.radians(bpnu))],[0,np.sin(np.radians(bpnu)),np.cos(np.radians(bpnu))]])
    M_bp=np.array([[np.cos(np.radians(bpmu)),0,np.sin(np.radians(bpmu))],[0,1,0],[-np.sin(np.radians(bpmu)),0,np.cos(np.radians(bpmu))]])
    UBt = 2 * np.pi * UB
    # Qv_pbの計算
    Qv_bp=UBt@(np.array([bp[0], bp[1], bp[2]]))
    Qv_bp = Qv_bp.reshape(-1, 1)  # 列ベクトルに変
    
    # fitting process
    def Omera_toration(omega):
        return np.array([[np.cos(np.radians(omega)),-np.sin(np.radians(omega)),0],[np.sin(np.radians(omega)),np.cos(np.radians(omega)),0],[0,0,1]])
    def N_rotation(nu):
        return np.array([[1,0,0],[0,np.cos(np.radians(nu)),-np.sin(np.radians(nu))],[0,np.sin(np.radians(nu)),np.cos(np.radians(nu))]])
    def M_rotation(mu):
        return np.array([[np.cos(np.radians(mu)),0,np.sin(np.radians(mu))],[0,1,0],[-np.sin(np.radians(mu)),0,np.cos(np.radians(mu))]])
    
    def objective0(angles):
        omega, mu, nu = angles
        rotation_matrix0 = Omera_toration(omega) @ N_rotation(mu) @ M_rotation(nu)
        transformed_vector0 = rotation_matrix0 @ Qv_bp
        return np.linalg.norm(transformed_vector0 - Qtheta_bp)

    # 最適化によって omega, mu, nu を求める
    initial_guess0 = [0, bpmu, bpnu]  # 初期値（全ての角度をゼロから開始）
    result0 = minimize(objective0, initial_guess0, method='BFGS')
    
    # 結果を表示
    omega0, mu0, nu0 = result0.x
    
    # s_iniの計算
    s_ini = omega0 + theta_bp
    # offsetの計算
    offset = bpc2 - s_ini
    
    # hklの位置から各angleを算出
    hkl_cal=cp[0]*astar+cp[1]*bstar+cp[2]*cstar
    Nhkl_cal=np.linalg.norm(hkl_cal)
    if fixe==0: # ei fix
        Ei = bpe
        Ef = bpe - cphw
    elif fixe==1: # ef fix
        Ei = bpe + cphw
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
    Qv_cal = UBt@(np.array([cp[0], cp[1], cp[2]]))
    Qv_cal=Qv_cal.reshape(-1, 1)
    
    # fitting process
    def objective(angles):
        omega, mu, nu = angles
        rotation_matrix = Omera_toration(omega) @ N_rotation(mu) @ M_rotation(nu)
        transformed_vector = rotation_matrix @ Qv_cal
        return np.linalg.norm(transformed_vector - Qtheta_cal)

    # 最適化によって omega, mu, nu を求める
    initial_guess = [0, 0, 0]  # 初期値（全ての角度をゼロから開始）
    result = minimize(objective, initial_guess, method='BFGS')
    
    # 結果を表示
    omega, mu, nu = result.x
    
    s_cal=omega+theta_cal
    omega_inst=s_cal+offset
    
    # アナライザとモノクロメータ
    d = 3.355  # PGの場合

    # C1とA1の計算
    C1 = np.degrees(np.arcsin((2 * np.pi / d) / (2 * np.sqrt(Ei / 2.072))))
    A1 = 2 * C1
    
    # C3とA3の計算
    C3 = np.degrees(np.arcsin((2 * np.pi / d) / (2 * np.sqrt(Ef / 2.072))))
    A3 = 2 * C3
        
    # 結果を辞書としてまとめる
    result = {
        'C1' : C1,
        'A1' : A1,
        'C2' : omega_inst,
        'A2' : phi_cal,
        'C3' : C3,
        'A3' : A3,
        'mu' : mu+bpmu,
        'nu' : nu+bpnu,
    }
    
    return result