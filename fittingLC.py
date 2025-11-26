import numpy as np
from scipy.optimize import minimize
import math

from RLcalc import RL_calc  # RL_calcを外部モジュールからインポート

# calculate_a2関数
def calculate_a2(ki, kf, N_hkl_bp):
    return np.degrees(np.arccos((ki**2 + kf**2 - N_hkl_bp**2) / (2 * ki * kf)))

# 自由度の判定
def determine_fit_params(crystal_system, hkl, initial_params=None):
    h, k, l = hkl
    # 初期値が渡されていない場合のデフォルト
    fixed_params = {
        "a": 5.0,
        "b": 5.0,
        "c": 5.0,
        "alpha": 90.0,
        "beta": 90.0,
        "gamma": 90.0,
    }
    
    # 初期値が指定されている場合は上書き
    if initial_params:
        for key in fixed_params:
            if key in initial_params:
                fixed_params[key] = initial_params[key]

    fit_params = []
    lattice_constraints = {}

    # 結晶系に応じた初期制約条件
    if crystal_system == "cubic":
        lattice_constraints = {"b": "a", "c": "a"}
        fixed_params["alpha"] = 90.0
        fixed_params["beta"] = 90.0
        fixed_params["gamma"] = 90.0
        fit_params = ["a"] if h or k or l else []
    elif crystal_system == "tetragonal":
        lattice_constraints = {"b": "a"}
        fixed_params["alpha"] = 90.0
        fixed_params["beta"] = 90.0
        fixed_params["gamma"] = 90.0
        fit_params = ["a"] if h or k else []
        if l:
            fit_params.append("c")
    elif crystal_system == "orthorhombic":
        fixed_params["alpha"] = 90.0
        fixed_params["beta"] = 90.0
        fixed_params["gamma"] = 90.0
        if h:
            fit_params.append("a")
        if k:
            fit_params.append("b")
        if l:
            fit_params.append("c")
    elif crystal_system == "hexagonal":
        fixed_params["alpha"] = 90.0
        fixed_params["beta"] = 90.0
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
def fit_lattice_constants(ki, kf, hkl, a2_measured, crystal_system, initial_params=None):
    # 初期値が与えられた場合、それを反映
    fixed_params, fit_params, lattice_constraints = determine_fit_params(crystal_system, hkl, initial_params)

    # フィッティング可能なパラメータがない場合
    if not fit_params:
        #print("No degrees of freedom available for fitting, returning initial values.")
        return fixed_params

    initial_guess = [fixed_params[param] for param in fit_params]

    # フィッティングを実行
    result = minimize(
        objective,
        initial_guess,
        args=(ki, kf, hkl, a2_measured, fit_params, fixed_params, lattice_constraints),
        method='Nelder-Mead'
    )
    
    if result.success:
        fitted_params = {param: value for param, value in zip(fit_params, result.x)}
        final_params = fixed_params.copy()
        final_params.update(fitted_params)

        # 格子定数の制約を適用
        for constraint, target in lattice_constraints.items():
            final_params[constraint] = final_params[target]

        return final_params
    else:
        raise ValueError("The fitting did not converge.")
