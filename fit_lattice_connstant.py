import numpy as np
from scipy.optimize import minimize
import math

from RLcalc import RL_calc  # RL_calcを外部モジュールからインポート
"""
# RL_calc関数
def RL_calc(a, b, c, alpha, beta, gamma):
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
    n_a = 2 * math.pi / V * b * c * math.sin(math.radians(alpha))
    n_b = 2 * math.pi / V * a * c * math.sin(math.radians(beta))
    n_c = 2 * math.pi / V * a * b * math.sin(math.radians(gamma))
    astar = n_a * np.array([1, 0, 0])
    bstar = n_b * np.array([
        math.cos(math.radians(gamma_star)),
        math.sin(math.radians(gamma_star)),
        0
    ])
    cstar = n_c * np.array([
        math.cos(math.radians(beta_star)),
        (math.cos(math.radians(alpha_star)) - math.cos(math.radians(beta_star)) * math.cos(math.radians(gamma_star))) / math.sin(math.radians(gamma_star)),
        V0 / math.sin(math.radians(gamma_star))
    ])
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
"""
# calculate_a2関数
def calculate_a2(ki, kf, N_hkl_bp):
    return np.degrees(np.arccos((ki**2 + kf**2 - N_hkl_bp**2) / (2 * ki * kf)))

# 自由度の判定
def determine_fit_params(crystal_system, hkl):
    h, k, l = hkl
    fixed_params = {
        "a": 5.0,
        "b": 5.0,
        "c": 5.0,
        "alpha": 90.0,
        "beta": 90.0,
        "gamma": 90.0,
    }
    fit_params = []
    lattice_constraints = {}

    # 結晶系に応じた初期制約条件
    if crystal_system == "cubic":
        lattice_constraints = {"b": "a", "c": "a"}
        fit_params = ["a"] if h or k or l else []
    elif crystal_system == "tetragonal":
        lattice_constraints = {"b": "a"}
        fit_params = ["a"] if h or k else []
        if l:
            fit_params.append("c")
    elif crystal_system == "orthorhombic":
        if h:
            fit_params.append("a")
        if k:
            fit_params.append("b")
        if l:
            fit_params.append("c")
    elif crystal_system == "hexagonal":
        fixed_params["gamma"] = 120.0
        lattice_constraints = {"b": "a"}
        fit_params = ["a"] if h or k else []
        if l:
            fit_params.append("c")
    elif crystal_system == "monoclinic":
        fixed_params["alpha"] = 90.0
        fixed_params["gamma"] = 90.0
        if h:
            fit_params.append("a")
        if k:
            fit_params.append("b")
        if l:
            fit_params.append("c")
        fit_params.append("beta")
    elif crystal_system == "triclinic":
        if h:
            fit_params.append("a")
        if k:
            fit_params.append("b")
        if l:
            fit_params.append("c")
        fit_params.extend(["alpha", "beta", "gamma"])

    return fixed_params, fit_params, lattice_constraints

# objective関数
def objective(params, ki, kf, hkl, a2_measured, fit_params, fixed_params, lattice_constraints):
    lattice_params = fixed_params.copy()
    for param, value in zip(fit_params, params):
        lattice_params[param] = value

    # 格子定数の制約を適用
    for constraint, target in lattice_constraints.items():
        lattice_params[constraint] = lattice_params[target]

    rl_values = RL_calc(
        lattice_params["a"], lattice_params["b"], lattice_params["c"],
        lattice_params["alpha"], lattice_params["beta"], lattice_params["gamma"]
    )
    astar = rl_values['astar']
    bstar = rl_values['bstar']
    cstar = rl_values['cstar']
    h, k, l = hkl
    hkl_bp = h * astar + k * bstar + l * cstar
    N_hkl_bp = np.linalg.norm(hkl_bp)

    a2_calc = calculate_a2(ki, kf, N_hkl_bp)
    return (a2_calc - a2_measured)**2

# fit_lattice_constants関数
def fit_lattice_constants(ki, kf, hkl, a2_measured, crystal_system):
    fixed_params, fit_params, lattice_constraints = determine_fit_params(crystal_system, hkl)

    # フィッティング可能なパラメータがない場合
    if not fit_params:
        print("フィッティング可能な自由度がありません。初期値を返します。")
        return fixed_params

    initial_guess = [fixed_params[param] for param in fit_params]

    result = minimize(objective, initial_guess, args=(ki, kf, hkl, a2_measured, fit_params, fixed_params, lattice_constraints), method='Nelder-Mead')
    if result.success:
        fitted_params = {param: value for param, value in zip(fit_params, result.x)}
        final_params = fixed_params.copy()
        final_params.update(fitted_params)

        # 格子定数の制約を適用
        for constraint, target in lattice_constraints.items():
            final_params[constraint] = final_params[target]

        return final_params
    else:
        raise ValueError("フィッティングが収束しませんでした。")

# 初期条件と実行
Ei=5
Ef=5
ki=(Ei/2.072)**(1/2)
kf=(Ef/2.072)**(1/2)
hkl = [1, 0, 0]
a2_measured = 50.0
crystal_system = "hexagonal"

fitted_params = fit_lattice_constants(ki, kf, hkl, a2_measured, crystal_system)
print(f"フィッティング結果: {fitted_params}")
