import math
import numpy as np
from scipy.interpolate import interp1d

def constqtime(flux_Ei, flux_cps, bpe, fixe, hw_ini, hw_fin, hw_inc, mcu):
    hw_tab = np.arange(hw_ini, hw_fin, hw_inc)
    if hw_tab[-1] != hw_fin:  # 最後の点が hw_fin でない場合
        hw_tab = np.append(hw_tab, hw_fin)
    
    if fixe==0: # ei fix
        Ef = np.zeros(len(hw_tab))
        for i in range(len(hw_tab)):
            Ei = bpe
            Ef[i] = bpe - hw_tab[i]
    elif fixe==1: # ef fix
        Ei = np.zeros(len(hw_tab))
        for i in range(len(hw_tab)):
            Ei[i] = bpe + hw_tab[i]
            Ef = bpe
    
    TE=np.zeros(len(hw_tab))
    
    if fixe==0: # ei fix
        # 線形補間関数を作成
        interp_func = interp1d(flux_Ei, flux_cps, kind='linear', fill_value='extrapolate')
        # Eiの補間値を計算
        cps_ei = interp_func(np.array([Ei]))
        cps_ef=np.zeros(len(hw_tab))
        for i in range(len(hw_tab)):
            # ターゲットXに対する補間値を計算
            cps_ef[i] = interp_func(Ef[i])
            TE[i] = mcu * cps_ei / cps_ef[i]
    elif fixe==1: # ef fix
        # 線形補間関数を作成
        interp_func = interp1d(flux_Ei, flux_cps, kind='linear', fill_value='extrapolate')
        # EiもしくはEfの補間値を計算
        cps_ef = interp_func(np.array([Ef]))
        cps_ei=np.zeros(len(hw_tab))
        for i in range(len(hw_tab)):
            # ターゲットXに対する補間値を計算
            cps_ei[i] = interp_func(Ei[i])
            TE[i] = mcu * cps_ef / cps_ei[i]
   
    sum_TE=np.sum(TE)
            
    # 60で割った商と余りを計算
    hour, min = divmod(sum_TE, 60)
   
    # 結果を辞書としてまとめる
    result = {
        'hour': hour,
        'min': min,
    }
    
    return result

def constetime(flux_Ei, flux_cps, bpe, fixe, hw_cal,h_ini,k_ini,l_ini,h_fin,k_fin,l_fin,h_inc,k_inc,l_inc, mcu):
    # 変数はhkl
    if h_ini != h_fin:
        h_tab = np.linspace(h_ini, h_fin, int((h_fin - h_ini) / h_inc) + 1)
    else:
        if k_ini != k_fin:
            h_tab = np.full(len(np.linspace(k_ini, k_fin, int((k_fin - k_ini) / k_inc) + 1)), h_ini)
        elif l_ini != l_fin:
            h_tab = np.full(len(np.linspace(l_ini, l_fin, int((l_fin - l_ini) / l_inc) + 1)), h_ini)

    if k_ini != k_fin:
        k_tab = np.linspace(k_ini, k_fin, int((k_fin - k_ini) / k_inc) + 1)
    else:
        if h_ini != h_fin:
            k_tab = np.full(len(h_tab), k_ini)
        elif l_ini != l_fin:
            k_tab = np.full(len(np.linspace(l_ini, l_fin, int((l_fin - l_ini) / l_inc) + 1)), k_ini)

    if l_ini != l_fin:
        l_tab = np.linspace(l_ini, l_fin, int((l_fin - l_ini) / l_inc) + 1)
    else:
        if h_ini != h_fin:
            l_tab = np.full(len(h_tab), l_ini)
        elif k_ini != k_fin:
            l_tab = np.full(len(k_tab), l_ini)

    if fixe==0: # ei fix
        Ei = bpe
        Ef = bpe - hw_cal
    elif fixe==1: # ef fix
        Ei = bpe + hw_cal
        Ef = bpe
    
    if fixe==0: # ei fix
        # 線形補間関数を作成
        interp_func = interp1d(flux_Ei, flux_cps, kind='linear', fill_value='extrapolate')
        # Eiの補間値を計算
        cps_ei = interp_func(np.array([Ei]))
        cps_ef = interp_func(np.array([Ef]))
        
        TE = mcu * cps_ei / cps_ef * len(l_tab)
        
    elif fixe==1: # ef fix
        # 線形補間関数を作成
        interp_func = interp1d(flux_Ei, flux_cps, kind='linear', fill_value='extrapolate')
        # EiもしくはEfの補間値を計算
        cps_ei = interp_func(np.array([Ei]))
        cps_ef = interp_func(np.array([Ef]))
        
        TE = mcu * cps_ef / cps_ei * len(l_tab)
    print(len(l_tab),TE)
    # 60で割った商と余りを計算
    hour, min = divmod(TE, 60)
   
    # 結果を辞書としてまとめる
    result = {
        'hour': hour,
        'min': min,
    }
    
    return result