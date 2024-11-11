import math
import numpy as np

def RL_calc(a, b, c, alpha, beta, gamma):
    """逆格子ベクトルとそのパラメータを計算し、辞書で返す"""
    # 体積 V と V0 の計算
    V = a * b * c * math.sqrt(
        1 - math.cos(math.radians(alpha))**2
        - math.cos(math.radians(beta))**2
        - math.cos(math.radians(gamma))**2
        + 2 * math.cos(math.radians(alpha)) 
        * math.cos(math.radians(beta)) 
        * math.cos(math.radians(gamma))
    )
    
    V0 = math.sqrt(
        1 - math.cos(math.radians(alpha))**2
        - math.cos(math.radians(beta))**2
        - math.cos(math.radians(gamma))**2
        + 2 * math.cos(math.radians(alpha)) 
        * math.cos(math.radians(beta)) 
        * math.cos(math.radians(gamma))
    )
    
    # 逆格子ベクトルの計算
    astar = a * np.array([1, 0, 0])
    bstar = b * np.array([
        math.cos(math.radians(gamma)),
        math.sin(math.radians(gamma)),
        0
    ])
    cstar = c * np.array([
        math.cos(math.radians(beta)),
        (math.cos(math.radians(alpha)) - math.cos(math.radians(beta)) * math.cos(math.radians(gamma))) / math.sin(math.radians(gamma)),
        V0 / math.sin(math.radians(gamma))
    ])
    
    # 角度の計算を簡潔化
    def acosd(x):
        return math.degrees(math.acos(x))

    alpha_star = acosd(
        (math.cos(math.radians(beta)) * math.cos(math.radians(gamma)) - math.cos(math.radians(alpha))) /
        (math.sin(math.radians(beta)) * math.sin(math.radians(gamma)))
    )
    
    beta_star = acosd(
        (math.cos(math.radians(alpha)) * math.cos(math.radians(gamma)) - math.cos(math.radians(beta))) /
        (math.sin(math.radians(alpha)) * math.sin(math.radians(gamma)))
    )
    
    gamma_star = acosd(
        (math.cos(math.radians(alpha)) * math.cos(math.radians(beta)) - math.cos(math.radians(gamma))) /
        (math.sin(math.radians(alpha)) * math.sin(math.radians(beta)))
    )
    
    # 逆格子定数の計算
    n_a = 2 * math.pi / V * b * c * math.sin(math.radians(alpha))
    n_b = 2 * math.pi / V * a * c * math.sin(math.radians(beta))
    n_c = 2 * math.pi / V * a * b * math.sin(math.radians(gamma))
    
    # 結果を辞書としてまとめる
    result = {
        'astar': astar,
        'bstar': bstar,
        'cstar': cstar,
        'alpha_star': alpha_star,
        'beta_star': beta_star,
        'gamma_star': gamma_star,
        'n_a': n_a,
        'n_b': n_b,
        'n_c': n_c,
        'V': V,
        'V0': V0
    }
    
    return result
