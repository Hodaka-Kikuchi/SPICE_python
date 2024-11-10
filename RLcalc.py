# 逆格子の計算。
import math
import numpy as np

def RL_calc(a,b,c,alpha,beta,gamma):
    V=a*b*c*math.sqrt(1-math.cos(math.radians(alpha))**2
                      -math.cos(math.radians(beta))**2
                      -math.cos(math.radians(gamma))**2
                      +2*math.cos(math.radians(alpha))*math.cos(math.radians(beta))*math.cos(math.radians(gamma)))
    
    
    V0=math.sqrt(1 - math.cos(math.radians(alpha))**2 
                   - math.cos(math.radians(beta))**2 
                   - math.cos(math.radians(gamma))**2 
                   + 2 * math.cos(math.radians(alpha)) 
                   * math.cos(math.radians(beta)) 
                   * math.cos(math.radians(gamma)))
    
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
    
    def acosd(x):
        return math.degrees(math.acos(x))
    
    alpha_star = acosd((math.cos(math.radians(beta)) * math.cos(math.radians(gamma)) - math.cos(math.radians(alpha))) /
               (math.sin(math.radians(beta)) * math.sin(math.radians(gamma))))
    
    beta_star = acosd((math.cos(math.radians(alpha)) * math.cos(math.radians(gamma)) - math.cos(math.radians(beta))) /
               (math.sin(math.radians(alpha)) * math.sin(math.radians(gamma))))
    
    gamma_star = acosd((math.cos(math.radians(alpha)) * math.cos(math.radians(beta)) - math.cos(math.radians(gamma))) /
               (math.sin(math.radians(alpha)) * math.sin(math.radians(beta))))
    
    n_a = 2 * math.pi / V * b * c * math.sin(math.radians(alpha))
    n_b = 2 * math.pi / V * a * c * math.sin(math.radians(beta))
    n_c = 2 * math.pi / V * a * b * math.sin(math.radians(gamma))

    n_V = n_a * n_b * n_c * V0
    
    # 結果を返す
    return astar,bstar,cstar,alpha_star,beta_star,gamma_star,n_a,n_b,n_c
    