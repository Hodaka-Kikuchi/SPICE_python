import math
import numpy as np
from scipy.optimize import minimize
import configparser
import os
import sys


def calculate_angle_with_plane(base_vec1, base_vec2, normal_vec, target_vec):
    """
    平面と任意のベクトルのなす角度を計算
    
    Parameters:
        base_vec1: 平面の1つ目の基準ベクトル
        base_vec2: 平面の2つ目の基準ベクトル
        normal_vec: 回転方向判定のための基準法線ベクトル
        target_vec: 計算対象のベクトル
    
    Returns:
        符号付き角度（度数単位）
    """
    # 平面の法線ベクトルを計算
    plane_normal = np.cross(base_vec1, base_vec2)
    plane_normal /= np.linalg.norm(plane_normal)  # 正規化
    
    # target_vecを平面に射影
    target_proj = target_vec - np.dot(target_vec, plane_normal) * plane_normal
    target_proj /= np.linalg.norm(target_proj)  # 射影ベクトルを正規化
    
    # base_vec1を基準に角度を計算
    dot_product = np.dot(base_vec1, target_proj)
    angle = -np.degrees(np.arccos(np.clip(dot_product, -1.0, 1.0)))
    
    # 外積で方向を判定
    cross_product = np.cross(base_vec1, target_proj)
    if np.dot(cross_product, normal_vec) < 0:
        angle = -angle  # 符号を付与
    
    return angle


def initial_value_with_multiple_planes(U, astar, bstar, cstar, hkl):
    """
    hklのベクトルと複数の平面のなす角度を計算
    
    Parameters:
        U: 回転行列
        astar, bstar, cstar: 基準ベクトル
        hkl: 計算対象のベクトル
    
    Returns:
        平面ごとの符号付き角度（度数単位）
    """
    vec1 = U[0, 0] * astar + U[0, 1] * bstar + U[0, 2] * cstar
    vec2 = U[1, 0] * astar + U[1, 1] * bstar + U[1, 2] * cstar
    vec3 = U[2, 0] * astar + U[2, 1] * bstar + U[2, 2] * cstar

    # 正規化
    Vec1 = vec1 / np.linalg.norm(vec1)
    Vec2 = vec2 / np.linalg.norm(vec2)
    Vec3 = vec3 / np.linalg.norm(vec3)

    # hklのベクトルを計算
    vect = hkl[0] * astar + hkl[1] * bstar + hkl[2] * cstar

    # 各平面との角度を計算
    angles_omega = calculate_angle_with_plane(Vec1, Vec2, Vec3, vect)
    angles_mu = calculate_angle_with_plane(Vec2, Vec3, Vec1, vect)
    angles_nu = -calculate_angle_with_plane(Vec1, Vec3, Vec2, vect)  # 符号が逆になるので補正

    # 角度を-90°～90°の範囲内に制限
    def normalize_angle(angle):
        if angle > 90:
            angle -= 180
        elif angle < -90:
            angle += 180
        return angle
    # angles_omegaがnanの場合、0に設定
    if np.isnan(angles_omega):
        angles_omega = 0
    angles_mu = normalize_angle(angles_mu)
    angles_nu = normalize_angle(angles_nu)

    return angles_omega, angles_mu, angles_nu

