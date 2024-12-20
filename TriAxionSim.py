#cd C:\DATA_HK\python\SPICE_python
# 右上にバージョン情報を表示
__version__ = '1.5.2'
"""
セマンティック バージョニング (Semantic Versioning)
セマンティック バージョニング（セムバ―、SemVer）は、バージョン番号を「MAJOR.MINOR.PATCH」の形式で表します。それぞれの部分には以下のような意味があります：

MAJOR（メジャー）バージョン:

後方互換性のない変更が導入された場合に増加します。
例：既存のAPIの変更、破壊的な変更。
MINOR（マイナー）バージョン

後方互換性のある新機能が追加された場合に増加します。
例：新しい機能の追加、既存機能の改良（後方互換性がある場合）。
PATCH（パッチ）バージョン:

後方互換性のあるバグ修正が行われた場合に増加します。
例：バグの修正、セキュリティ修正。

バージョン番号の例
1.0.0: 最初の安定版リリース。
1.1.0: 後方互換性のある新機能が追加されたリリース。
1.1.1: バグ修正やマイナーな改良が行われたリリース。
2.0.0: 後方互換性のない変更が導入されたリリース。

新機能付与とバグ修正を行った場合、patchバージョンを+、ラストの数値は0に戻す(更新の必要はない。)。
常に左側の数値の更新が優先。
"""

# tlinterのインポート
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog

# 確定的progressbarの設置
from tkinter import messagebox
from tkinter import simpledialog
import pyperclip

import configparser
# osのインポート
import os

# ギリシャ文字の定義
import sympy as sm
sm.init_printing()

mu    = sm.Symbol("μ")
nu    = sm.Symbol("ν")
theta = sm.Symbol("θ")# "\theta"から変更
alpha    = sm.Symbol("α")
beta    = sm.Symbol("β")
gamma    = sm.Symbol("γ")
AA    = sm.Symbol("Å")

# 数値計算を行うためのパッケージのインポート
from scipy.interpolate import interp1d

# smoothing処理のパッケージ
from scipy.ndimage import gaussian_filter

import numpy as np
import math
import matplotlib
matplotlib.use('TkAgg')#よくわかんないけどこれないとexe化したときにグラフが表示されない。超重要
import matplotlib.pyplot as plt
import matplotlib.widgets as wg #from  matplotlib.widgets import Slider
from mpl_toolkits.mplot3d import Axes3D

from itertools import product

import matplotlib.ticker as mticker

from matplotlib.backend_tools import ToolBase, ToolToggleBase

from matplotlib.colors import LogNorm

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from matplotlib.widgets import TextBox

from matplotlib.patches import Rectangle
import matplotlib.patches as patches

from fractions import Fraction

import statistics as stat
from scipy import stats
from scipy.stats import norm
from scipy.optimize import curve_fit

# ファイル読み込みのためのやつ
import re
import sys

# webに飛ぶやつ
import webbrowser

import csv
from tkinter import filedialog

from tkinter import PhotoImage

#windowの作成
root=tk.Tk()
#windowのタイトル変更
root.title(f"TriAxionSim ver: {__version__}")
# TriAxionSim: 三軸 (triple-axis) と「軌跡」や「軸」 (axion) 、Simulationを意識
#windowのサイズ指定
root.geometry("550x840")#550*840

# ロゴを設定
# 実行時のリソースパスを設定
def resource_path(relative_path):
    """PyInstallerでパスを解決する関数"""
    if hasattr(sys, '_MEIPASS'):
        return os.path.join(sys._MEIPASS, relative_path)
    return os.path.join(os.path.abspath("."), relative_path)

# アイコンの設定
logo_path = resource_path("logo2.ico")
root.iconbitmap(logo_path)

# フォント設定
default_font_size = 10
entry_font = ('Helvetica', default_font_size)
#ttk.entry(フレーム名,text=*****,, font=entry_font)で統一可能。ラベルも同様。

# ×ボタンを押すと作成したグラフ全てがクリアされる。
# ウィンドウが閉じられたときの処理
def on_closing():
    # すべてのグラフウィンドウを閉じる
    plt.close('all')
    root.destroy()

root.protocol("WM_DELETE_WINDOW", on_closing)  # ウィンドウが閉じられるときの振る舞いを指定

# iniファイルの読み込み
def load_values_from_ini():
    config = configparser.ConfigParser()
    # .exe化した場合に対応する
    if getattr(sys, 'frozen', False):
        # .exeの場合、sys.argv[0]が実行ファイルのパスになる
        ini_path = os.path.join(os.path.dirname(sys.argv[0]), 'config.ini')
    else:
        # .pyの場合、__file__がスクリプトのパスになる
        ini_path = os.path.join(os.path.dirname(__file__), 'config.ini')
    config.read(ini_path)
    
    # 各エントリに対応する値を読み込み、挿入
    la.delete(0, tk.END)  # 既存の値をクリア
    la.insert(0, config['sample'].get('a', '5.026307'))
    lb.delete(0, tk.END)  # 既存の値をクリア
    lb.insert(0, config['sample'].get('b', '5.026307'))
    lc.delete(0, tk.END)  # 既存の値をクリア
    lc.insert(0, config['sample'].get('c', '13.784500'))
    lc_alpha.delete(0, tk.END)  # 既存の値をクリア
    lc_alpha.insert(0, config['sample'].get('alpha', '90'))
    lc_beta.delete(0, tk.END)  # 既存の値をクリア
    lc_beta.insert(0, config['sample'].get('beta', '90'))
    lc_gamma.delete(0, tk.END)  # 既存の値をクリア
    lc_gamma.insert(0, config['sample'].get('gamma', '120'))
    
    sv1_h.delete(0, tk.END)  # 既存の値をクリア
    sv1_h.insert(0, config['sample'].get('h1', '1'))
    sv1_k.delete(0, tk.END)  # 既存の値をクリア
    sv1_k.insert(0, config['sample'].get('k1', '0'))
    sv1_l.delete(0, tk.END)  # 既存の値をクリア
    sv1_l.insert(0, config['sample'].get('l1', '0'))
    sv2_h.delete(0, tk.END)  # 既存の値をクリア
    sv2_h.insert(0, config['sample'].get('h2', '0'))
    sv2_k.delete(0, tk.END)  # 既存の値をクリア
    sv2_k.insert(0, config['sample'].get('k2', '1'))
    sv2_l.delete(0, tk.END)  # 既存の値をクリア
    sv2_l.insert(0, config['sample'].get('l2', '0'))
    
    Energy.delete(0, tk.END)  # 既存の値をクリア
    Energy.insert(0, config['BPP'].get('Energy', '4.8'))
    bp_h.delete(0, tk.END)  # 既存の値をクリア
    bp_h.insert(0, config['BPP'].get('bp_h', '1.5'))
    bp_k.delete(0, tk.END)  # 既存の値をクリア
    bp_k.insert(0, config['BPP'].get('bp_k', '0'))
    bp_l.delete(0, tk.END)  # 既存の値をクリア
    bp_l.insert(0, config['BPP'].get('bp_l', '0'))
    bp_c2.delete(0, tk.END)  # 既存の値をクリア
    bp_c2.insert(0, config['BPP'].get('bp_c2', '49.75'))
    bp_mu.delete(0, tk.END)  # 既存の値をクリア
    bp_mu.insert(0, config['BPP'].get('bp_mu', '0'))
    bp_nu.delete(0, tk.END)  # 既存の値をクリア
    bp_nu.insert(0, config['BPP'].get('bp_nu', '0'))
    
    # ラジオボタンの初期状態を読み込む
    eief_value = int(config['instrument'].get('eief', '1'))  # 1がデフォルト
    # eiefの初期値を設定
    eief.set(eief_value)
    
    # ラベルの自動更新
    update_label()
    # collimator & mosaic
    div_1st_m.delete(0, tk.END)  # 既存の値をクリア
    div_1st_m.insert(0, config['instrument'].get('div_1st_m', '1.2'))
    div_2nd_h.delete(0, tk.END)  # 既存の値をクリア
    div_2nd_h.insert(0, config['instrument'].get('div_2nd_h', '120'))
    div_2nd_v.delete(0, tk.END)  # 既存の値をクリア
    div_2nd_v.insert(0, config['instrument'].get('div_2nd_v', '200'))
    div_3rd_h.delete(0, tk.END)  # 既存の値をクリア
    div_3rd_h.insert(0, config['instrument'].get('div_3rd_h', '80'))
    div_3rd_v.delete(0, tk.END)  # 既存の値をクリア
    div_3rd_v.insert(0, config['instrument'].get('div_3rd_v', '200'))
    div_4th_h.delete(0, tk.END)  # 既存の値をクリア
    div_4th_h.insert(0, config['instrument'].get('div_4th_h', '120'))
    div_4th_v.delete(0, tk.END)  # 既存の値をクリア
    div_4th_v.insert(0, config['instrument'].get('div_4th_v', '200'))
    mos_sam_h.delete(0, tk.END)  # 既存の値をクリア
    mos_sam_h.insert(0, config['sample'].get('mos_sam_h', '50'))
    mos_sam_v.delete(0, tk.END)  # 既存の値をクリア
    mos_sam_v.insert(0, config['sample'].get('mos_sam_v', '50'))
    
    mos_mono_h.delete(0, tk.END)  # 既存の値をクリア
    mos_mono_h.insert(0, config['instrument'].get('mos_mono_h', '50'))
    mos_mono_v.delete(0, tk.END)  # 既存の値をクリア
    mos_mono_v.insert(0, config['instrument'].get('mos_mono_v', '50'))
    mos_ana_h.delete(0, tk.END)  # 既存の値をクリア
    mos_ana_h.insert(0, config['instrument'].get('mos_ana_h', '50'))
    mos_ana_v.delete(0, tk.END)  # 既存の値をクリア
    mos_ana_v.insert(0, config['instrument'].get('mos_ana_v', '50'))
    
    # hardware limit
    hwl2f.delete(0, tk.END)  # 既存の値をクリア
    hwl2f.insert(0, config['instrument'].get('maxC1', '19.9305'))
    hwl2t.delete(0, tk.END)  # 既存の値をクリア
    hwl2t.insert(0, config['instrument'].get('minC1', '58.482'))
    hwl3f.delete(0, tk.END)  # 既存の値をクリア
    hwl3f.insert(0, config['instrument'].get('maxA1', '19.9305'))
    hwl3t.delete(0, tk.END)  # 既存の値をクリア
    hwl3t.insert(0, config['instrument'].get('minA1', '58.482'))
    hwl4f.delete(0, tk.END)  # 既存の値をクリア
    hwl4f.insert(0, config['instrument'].get('maxC2', '-180'))
    hwl4t.delete(0, tk.END)  # 既存の値をクリア
    hwl4t.insert(0, config['instrument'].get('minC2', '180'))
    hwl5f.delete(0, tk.END)  # 既存の値をクリア
    hwl5f.insert(0, config['instrument'].get('maxA2', '6'))
    hwl5t.delete(0, tk.END)  # 既存の値をクリア
    hwl5t.insert(0, config['instrument'].get('minA2', '120'))
    hwl6f.delete(0, tk.END)  # 既存の値をクリア
    hwl6f.insert(0, config['instrument'].get('maxC3', '19.9305'))
    hwl6t.delete(0, tk.END)  # 既存の値をクリア
    hwl6t.insert(0, config['instrument'].get('minC3', '58.482'))
    hwl7f.delete(0, tk.END)  # 既存の値をクリア
    hwl7f.insert(0, config['instrument'].get('maxA3', '39.861'))
    hwl7t.delete(0, tk.END)  # 既存の値をクリア
    hwl7t.insert(0, config['instrument'].get('minA3', '116.964'))
    hwl8f.delete(0, tk.END)  # 既存の値をクリア
    hwl8f.insert(0, config['instrument'].get('maxmu', '-5'))
    hwl8t.delete(0, tk.END)  # 既存の値をクリア
    hwl8t.insert(0, config['instrument'].get('minmu', '5'))
    hwl9f.delete(0, tk.END)  # 既存の値をクリア
    hwl9f.insert(0, config['instrument'].get('maxnu', '-5'))
    hwl9t.delete(0, tk.END)  # 既存の値をクリア
    hwl9t.insert(0, config['instrument'].get('minnu', '5'))
    acna.delete(0, tk.END)  # 既存の値をクリア
    acna.insert(0, config['instrument'].get('blade_num', '7'))
    
    #figure option
    # チェックボックス
    fig_spec_value = int(config['option'].get('fig_spec', '1'))  # 1がデフォルト
    # fig_specの初期値を設定
    fig_spec.set(fig_spec_value)
    fig_reso_value = int(config['option'].get('fig_reso', '1'))  # 1がデフォルト  
    # fig_resoの初期値を設定
    fig_reso.set(fig_reso_value)
    calc_hf_value = int(config['option'].get('calc_hf', '1'))  # 1がデフォルト
    # calc_hfの初期値を設定
    calc_hf.set(calc_hf_value)
    
    # fig_resoの初期値を設定
    fig_reci_value = int(config['option'].get('fig_reci', '1'))  # 1がデフォルト
    fig_reci.set(fig_reci_value)

def save_values_to_ini():
    """
    現在のウィジェットの値をINIファイルに保存する
    """
    config = configparser.ConfigParser()
    
    # INIファイルのパスを決定
    if getattr(sys, 'frozen', False):
        ini_path = os.path.join(os.path.dirname(sys.argv[0]), 'config.ini')
    else:
        ini_path = os.path.join(os.path.dirname(__file__), 'config.ini')
    
    # 既存のINIファイルを読み込む
    if os.path.exists(ini_path):
        config.read(ini_path, encoding='utf-8')  # UTF-8で読み込み
    
    # 必要なセクションがなければ作成
    if 'sample' not in config:
        config['sample'] = {}
    if 'instrument' not in config:
        config['instrument'] = {}

    # 'sample'セクションを更新
    config['sample'].update({
        'a': la.get(),
        'b': lb.get(),
        'c': lc.get(),
        'alpha': lc_alpha.get(),
        'beta': lc_beta.get(),
        'gamma': lc_gamma.get(),
        'h1': sv1_h.get(),
        'k1': sv1_k.get(),
        'l1': sv1_l.get(),
        'h2': sv2_h.get(),
        'k2': sv2_k.get(),
        'l2': sv2_l.get(),
        'mos_sam_h': mos_sam_h.get(),
        'mos_sam_v': mos_sam_v.get(),
    })
    
    # 'sample'セクションを更新
    config['BPP'].update({
        'Energy': Energy.get(),
        'bp_h': bp_h.get(),
        'bp_k': bp_k.get(),
        'bp_l': bp_l.get(),
        'bp_c2': bp_c2.get(),
        'bp_mu': bp_mu.get(),
        'bp_nu': bp_nu.get(),
    })

    # 'instrument'セクションを更新
    config['instrument'].update({
        'eief': str(eief.get()),  # ラジオボタンの状態を保存
        'div_1st_m': div_1st_m.get(),
        'div_2nd_h': div_2nd_h.get(),
        'div_2nd_v': div_2nd_v.get(),
        'div_3rd_h': div_3rd_h.get(),
        'div_3rd_v': div_3rd_v.get(),
        'div_4th_h': div_4th_h.get(),
        'div_4th_v': div_4th_v.get(),
        'mos_mono_h': mos_mono_h.get(),
        'mos_mono_v': mos_mono_v.get(),
        'mos_ana_h': mos_ana_h.get(),
        'mos_ana_v': mos_ana_v.get(),
        'maxC1': hwl2f.get(),
        'minC1': hwl2t.get(),
        'maxA1': hwl3f.get(),
        'minA1': hwl3t.get(),
        'maxC2': hwl4f.get(),
        'minC2': hwl4t.get(),
        'maxA2': hwl5f.get(),
        'minA2': hwl5t.get(),
        'maxC3': hwl6f.get(),
        'minC3': hwl6t.get(),
        'maxA3': hwl7f.get(),
        'minA3': hwl7t.get(),
        'maxmu': hwl8f.get(),
        'minmu': hwl8t.get(),
        'maxnu': hwl9f.get(),
        'minnu': hwl9t.get(),
        'blade_num': acna.get(),
    })
    
    # 'sample'セクションを更新
    config['option'].update({
        'fig_spec': str(fig_spec.get()),
        'fig_reso': str(fig_reso.get()),
        'fig_reci': str(fig_reci.get()),
        'calc_hf': str(calc_hf.get()),
    })
    
    # INIファイルに書き込み
    with open(ini_path, 'w') as configfile:
        config.write(configfile)

# fuction setting
from RLcalc import RL_calc  #
from UBcalc import UB_calc  #
from anglecalc import angle_calc    #
from anglecalc2 import angle_calc2    #
from anglecalc3 import angle_calc3    #
from specfigscan import plot_spectrometer #
from fittingLC import fit_lattice_constants
from QEresolution import calcresolution
from QEresolution_scan import calcresolution_scan
from fig_reciprocal_space import plot_reciprocal_space

from specfigscan_gif import plot_spectrometer_with_gif #画像をgif保存したいとき

# GUIの配分を決める。
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=2)
root.rowconfigure(1, weight=2)
root.rowconfigure(2, weight=4)
root.rowconfigure(3, weight=2)
root.rowconfigure(4, weight=5)
root.rowconfigure(5, weight=3)
root.rowconfigure(6, weight=3)

# ファイル選択のフレームの作成と設置
frame1 = ttk.Labelframe(root,text= "lattice information")
frame1.grid(row=0,column=0,sticky="NSEW")

frame1.columnconfigure(0, weight=1)
frame1.columnconfigure(1, weight=1)
frame1.columnconfigure(2, weight=1)
frame1.columnconfigure(3, weight=1)
frame1.columnconfigure(4, weight=1)
frame1.columnconfigure(5, weight=1)
frame1.rowconfigure(0, weight=1)
frame1.rowconfigure(1, weight=1)

def get_parameters():
    """GUI入力から格子パラメータを取得し辞書に格納"""
    parameters = {
        'a': float(la.get()),
        'b': float(lb.get()),
        'c': float(lc.get()),
        'alpha': float(lc_alpha.get()),
        'beta': float(lc_beta.get()),
        'gamma': float(lc_gamma.get())
    }
    return parameters

def on_Rlcalc():
    """RLtableを計算して辞書で返す"""
    params = get_parameters()
    RLtable = RL_calc(**params)  # RL_calcに辞書を展開して渡す
    return RLtable

def on_UBcalc():
    """UBtableを計算して返す"""
    params = get_parameters()
    
    # 散乱面の取得
    sv1 = np.array([float(sv1_h.get()), float(sv1_k.get()), float(sv1_l.get())])
    sv2 = np.array([float(sv2_h.get()), float(sv2_k.get()), float(sv2_l.get())])
    
    # RLtableを取得し、辞書から必要な変数を取り出す
    RLtable = on_Rlcalc()
    astar = RLtable['astar']
    bstar = RLtable['bstar']
    cstar = RLtable['cstar']
    alpha_star = RLtable['alpha_star']
    beta_star = RLtable['beta_star']
    gamma_star = RLtable['gamma_star']
    n_a = RLtable['n_a']
    n_b = RLtable['n_b']
    n_c = RLtable['n_c']
    
    #Sense=sense.get()
    
    # UBtableを計算
    UBtable = UB_calc(
        sv1, sv2, astar, bstar, cstar, alpha_star, beta_star, gamma_star, 
        n_a, n_b, n_c, **params
    )
    return UBtable

# 格子定数を入力する欄
lc1 = tk.Label(frame1,text='a (Å)')
lc1.grid(row=0, column=0,sticky="NSEW")
la = ttk.Entry(frame1)
la.grid(row=1, column=0,sticky="NSEW")
#la.insert(0,'5.026307')

lc2 = tk.Label(frame1,text='b (Å)')
lc2.grid(row=0, column=1,sticky="NSEW")
lb = ttk.Entry(frame1)
lb.grid(row=1, column=1,sticky="NSEW")
#lb.insert(0,'5.026307')

lc3 = tk.Label(frame1,text='c (Å)')
lc3.grid(row=0, column=2,sticky="NSEW")
lc = ttk.Entry(frame1)
lc.grid(row=1, column=2,sticky="NSEW")
#lc.insert(0,'13.784500')

lc4 = tk.Label(frame1,text='α (deg)')
lc4.grid(row=0, column=3,sticky="NSEW")
lc_alpha = ttk.Entry(frame1)
lc_alpha.grid(row=1, column=3,sticky="NSEW")
#lc_alpha.insert(0,'90')

lc5 = tk.Label(frame1,text='β (deg)')
lc5.grid(row=0, column=4,sticky="NSEW")
lc_beta = ttk.Entry(frame1)
lc_beta.grid(row=1, column=4,sticky="NSEW")
#lc_beta.insert(0,'90')

lc6 = tk.Label(frame1,text='γ (deg)')
lc6.grid(row=0, column=5,sticky="NSEW")
lc_gamma = ttk.Entry(frame1)
lc_gamma.grid(row=1, column=5,sticky="NSEW")
#lc_gamma.insert(0,'120')

# ファイル選択のフレームの作成と設置
frame2 = ttk.Labelframe(root,text= "scattering plane")
frame2.grid(row=1,column=0,sticky="NSEW")

frame2.columnconfigure(0, weight=1)
frame2.columnconfigure(1, weight=1)
frame2.rowconfigure(0, weight=1)

frame2a = ttk.Labelframe(frame2,text= "axis 1")
frame2a.grid(row=0,column=0,sticky="NSEW")

frame2a.columnconfigure(0, weight=1)
frame2a.columnconfigure(1, weight=1)
frame2a.columnconfigure(2, weight=1)
frame2a.rowconfigure(0, weight=1)
frame2a.rowconfigure(1, weight=1)

frame2b = ttk.Labelframe(frame2,text= "axis 2")
frame2b.grid(row=0,column=1,sticky="NSEW")

frame2b.columnconfigure(0, weight=1)
frame2b.columnconfigure(1, weight=1)
frame2b.columnconfigure(2, weight=1)
frame2b.rowconfigure(0, weight=1)
frame2b.rowconfigure(1, weight=1)

# 散乱面を入力する欄
sv1 = tk.Label(frame2a,text='h')
sv1.grid(row=0, column=0,sticky="NSEW")
sv1_h = ttk.Entry(frame2a)
sv1_h.grid(row=1, column=0,sticky="NSEW")
#sv1_h.insert(0,'1')

sv1 = tk.Label(frame2a,text='k')
sv1.grid(row=0, column=1,sticky="NSEW")
sv1_k = ttk.Entry(frame2a)
sv1_k.grid(row=1, column=1,sticky="NSEW")
#sv1_k.insert(0,'0')

sv1 = tk.Label(frame2a,text='l')
sv1.grid(row=0, column=2,sticky="NSEW")
sv1_l = ttk.Entry(frame2a)
sv1_l.grid(row=1, column=2,sticky="NSEW")
#sv1_l.insert(0,'0')

sv2 = tk.Label(frame2b,text='h')
sv2.grid(row=0, column=0,sticky="NSEW")
sv2_h = ttk.Entry(frame2b)
sv2_h.grid(row=1, column=0,sticky="NSEW")
#sv2_h.insert(0,'0')

sv2 = tk.Label(frame2b,text='k')
sv2.grid(row=0, column=1,sticky="NSEW")
sv2_k = ttk.Entry(frame2b)
sv2_k.grid(row=1, column=1,sticky="NSEW")
#sv2_k.insert(0,'1')

sv2 = tk.Label(frame2b,text='l')
sv2.grid(row=0, column=2,sticky="NSEW")
sv2_l = ttk.Entry(frame2b)
sv2_l.grid(row=1, column=2,sticky="NSEW")
#sv2_l.insert(0,'0')

# UB matrixの表示
frame3 = ttk.Labelframe(root,text= "matrix display")
frame3.grid(row=2,column=0,sticky="NSEW")

frame3.columnconfigure(0, weight=1)
frame3.columnconfigure(1, weight=1)
frame3.columnconfigure(2, weight=1)
frame3.rowconfigure(0, weight=3)
frame3.rowconfigure(1, weight=1)

frame3a = ttk.Labelframe(frame3,text= "U matrix")
frame3a.grid(row=0,column=0,sticky="NSEW")

frame3a.columnconfigure(0, weight=1)
frame3a.columnconfigure(1, weight=1)
frame3a.columnconfigure(2, weight=1)
frame3a.rowconfigure(0, weight=1)
frame3a.rowconfigure(1, weight=1)
frame3a.rowconfigure(2, weight=1)

frame3b = ttk.Labelframe(frame3,text= "B matrix")
frame3b.grid(row=0,column=1,sticky="NSEW")

frame3b.columnconfigure(0, weight=1)
frame3b.columnconfigure(1, weight=1)
frame3b.columnconfigure(2, weight=1)
frame3b.rowconfigure(0, weight=1)
frame3b.rowconfigure(1, weight=1)
frame3b.rowconfigure(2, weight=1)

frame3c = ttk.Labelframe(frame3,text= "UB matrix")
frame3c.grid(row=0,column=2,sticky="NSEW")

frame3c.columnconfigure(0, weight=1)
frame3c.columnconfigure(1, weight=1)
frame3c.columnconfigure(2, weight=1)
frame3c.rowconfigure(0, weight=1)
frame3c.rowconfigure(1, weight=1)
frame3c.rowconfigure(2, weight=1)

# Umatrixの表示
txt_u_11 = ttk.Entry(frame3a,state="readonly")
txt_u_11.insert(0,'0')
txt_u_11.grid(row=0, column=0,sticky="NSEW")
txt_u_12 = ttk.Entry(frame3a,state="readonly")
txt_u_12.insert(0,'0')
txt_u_12.grid(row=0, column=1,sticky="NSEW")
txt_u_13 = ttk.Entry(frame3a,state="readonly")
txt_u_13.insert(0,'0')
txt_u_13.grid(row=0, column=2,sticky="NSEW")
txt_u_21 = ttk.Entry(frame3a,state="readonly")
txt_u_21.insert(0,'0')
txt_u_21.grid(row=1, column=0,sticky="NSEW")
txt_u_22 = ttk.Entry(frame3a,state="readonly")
txt_u_22.insert(0,'0')
txt_u_22.grid(row=1, column=1,sticky="NSEW")
txt_u_23 = ttk.Entry(frame3a,state="readonly")
txt_u_23.insert(0,'0')
txt_u_23.grid(row=1, column=2,sticky="NSEW")
txt_u_31 = ttk.Entry(frame3a,state="readonly")
txt_u_31.insert(0,'0')
txt_u_31.grid(row=2, column=0,sticky="NSEW")
txt_u_32 = ttk.Entry(frame3a,state="readonly")
txt_u_32.insert(0,'0')
txt_u_32.grid(row=2, column=1,sticky="NSEW")
txt_u_33 = ttk.Entry(frame3a,state="readonly")
txt_u_33.insert(0,'0')
txt_u_33.grid(row=2, column=2,sticky="NSEW")

# Bmatrixの表示
txt_b_11 = ttk.Entry(frame3b,state="readonly")
txt_b_11.insert(0,'0')
txt_b_11.grid(row=0, column=0,sticky="NSEW")
txt_b_12 = ttk.Entry(frame3b,state="readonly")
txt_b_12.insert(0,'0')
txt_b_12.grid(row=0, column=1,sticky="NSEW")
txt_b_13 = ttk.Entry(frame3b,state="readonly")
txt_b_13.insert(0,'0')
txt_b_13.grid(row=0, column=2,sticky="NSEW")
txt_b_21 = ttk.Entry(frame3b,state="readonly")
txt_b_21.insert(0,'0')
txt_b_21.grid(row=1, column=0,sticky="NSEW")
txt_b_22 = ttk.Entry(frame3b,state="readonly")
txt_b_22.insert(0,'0')
txt_b_22.grid(row=1, column=1,sticky="NSEW")
txt_b_23 = ttk.Entry(frame3b,state="readonly")
txt_b_23.insert(0,'0')
txt_b_23.grid(row=1, column=2,sticky="NSEW")
txt_b_31 = ttk.Entry(frame3b,state="readonly")
txt_b_31.insert(0,'0')
txt_b_31.grid(row=2, column=0,sticky="NSEW")
txt_b_32 = ttk.Entry(frame3b,state="readonly")
txt_b_32.insert(0,'0')
txt_b_32.grid(row=2, column=1,sticky="NSEW")
txt_b_33 = ttk.Entry(frame3b,state="readonly")
txt_b_33.insert(0,'0')
txt_b_33.grid(row=2, column=2,sticky="NSEW")

# UBmatrixの表示
txt_ub_11 = ttk.Entry(frame3c,state="readonly")
txt_ub_11.insert(0,'0')
txt_ub_11.grid(row=0, column=0,sticky="NSEW")
txt_ub_12 = ttk.Entry(frame3c,state="readonly")
txt_ub_12.insert(0,'0')
txt_ub_12.grid(row=0, column=1,sticky="NSEW")
txt_ub_13 = ttk.Entry(frame3c,state="readonly")
txt_ub_13.insert(0,'0')
txt_ub_13.grid(row=0, column=2,sticky="NSEW")
txt_ub_21 = ttk.Entry(frame3c,state="readonly")
txt_ub_21.insert(0,'0')
txt_ub_21.grid(row=1, column=0,sticky="NSEW")
txt_ub_22 = ttk.Entry(frame3c,state="readonly")
txt_ub_22.insert(0,'0')
txt_ub_22.grid(row=1, column=1,sticky="NSEW")
txt_ub_23 = ttk.Entry(frame3c,state="readonly")
txt_ub_23.insert(0,'0')
txt_ub_23.grid(row=1, column=2,sticky="NSEW")
txt_ub_31 = ttk.Entry(frame3c,state="readonly")
txt_ub_31.insert(0,'0')
txt_ub_31.grid(row=2, column=0,sticky="NSEW")
txt_ub_32 = ttk.Entry(frame3c,state="readonly")
txt_ub_32.insert(0,'0')
txt_ub_32.grid(row=2, column=1,sticky="NSEW")
txt_ub_33 = ttk.Entry(frame3c,state="readonly")
txt_ub_33.insert(0,'0')
txt_ub_33.grid(row=2, column=2,sticky="NSEW")

# 1つのボタンで両方の計算を実行
def displayUB():
    # on_Rlcalc の結果を取得
    RLtable = on_Rlcalc()
    #print("RLtable:", RLtable)
    
    # on_UBcalc の結果を取得
    UBtable = on_UBcalc()
    #print("UBtable:", UBtable)
    
    # U matrixを表示
    txt_u_11.config(state="normal")  # 一時的に編集可能に
    txt_u_11.delete(0, tk.END)
    txt_u_11.insert(0, round(UBtable['U'][0,0],4))
    txt_u_11.config(state="readonly") # 編集不可に設定
    txt_u_12.config(state="normal")  # 一時的に編集可能に
    txt_u_12.delete(0, tk.END)
    txt_u_12.insert(0, round(UBtable['U'][0,1],4))
    txt_u_12.config(state="readonly") # 編集不可に設定
    txt_u_13.config(state="normal")  # 一時的に編集可能に
    txt_u_13.delete(0, tk.END)
    txt_u_13.insert(0, round(UBtable['U'][0,2],4))
    txt_u_13.config(state="readonly") # 編集不可に設定
    txt_u_21.config(state="normal")  # 一時的に編集可能に
    txt_u_21.delete(0, tk.END)
    txt_u_21.insert(0, round(UBtable['U'][1,0],4))
    txt_u_21.config(state="readonly") # 編集不可に設定
    txt_u_22.config(state="normal")  # 一時的に編集可能に
    txt_u_22.delete(0, tk.END)
    txt_u_22.insert(0, round(UBtable['U'][1,1],4))
    txt_u_22.config(state="readonly") # 編集不可に設定
    txt_u_23.config(state="normal")  # 一時的に編集可能に
    txt_u_23.delete(0, tk.END)
    txt_u_23.insert(0, round(UBtable['U'][1,2],4))
    txt_u_23.config(state="readonly") # 編集不可に設定
    txt_u_31.config(state="normal")  # 一時的に編集可能に
    txt_u_31.delete(0, tk.END)
    txt_u_31.insert(0, round(UBtable['U'][2,0],4))
    txt_u_31.config(state="readonly") # 編集不可に設定
    txt_u_32.config(state="normal")  # 一時的に編集可能に
    txt_u_32.delete(0, tk.END)
    txt_u_32.insert(0, round(UBtable['U'][2,1],4))
    txt_u_32.config(state="readonly") # 編集不可に設定
    txt_u_33.config(state="normal")  # 一時的に編集可能に
    txt_u_33.delete(0, tk.END)
    txt_u_33.insert(0, round(UBtable['U'][2,2],4))
    txt_u_33.config(state="readonly") # 編集不可に設定
    
    # B matrixを表示
    txt_b_11.config(state="normal")  # 一時的に編集可能に
    txt_b_11.delete(0, tk.END)
    txt_b_11.insert(0, round(UBtable['B'][0,0],4))
    txt_b_11.config(state="readonly") # 編集不可に設定
    txt_b_12.config(state="normal")  # 一時的に編集可能に
    txt_b_12.delete(0, tk.END)
    txt_b_12.insert(0, round(UBtable['B'][0,1],4))
    txt_b_12.config(state="readonly") # 編集不可に設定
    txt_b_13.config(state="normal")  # 一時的に編集可能に
    txt_b_13.delete(0, tk.END)
    txt_b_13.insert(0, round(UBtable['B'][0,2],4))
    txt_b_13.config(state="readonly") # 編集不可に設定
    txt_b_21.config(state="normal")  # 一時的に編集可能に
    txt_b_21.delete(0, tk.END)
    txt_b_21.insert(0, round(UBtable['B'][1,0],4))
    txt_b_21.config(state="readonly") # 編集不可に設定
    txt_b_22.config(state="normal")  # 一時的に編集可能に
    txt_b_22.delete(0, tk.END)
    txt_b_22.insert(0, round(UBtable['B'][1,1],4))
    txt_b_22.config(state="readonly") # 編集不可に設定
    txt_b_23.config(state="normal")  # 一時的に編集可能に
    txt_b_23.delete(0, tk.END)
    txt_b_23.insert(0, round(UBtable['B'][1,2],4))
    txt_b_23.config(state="readonly") # 編集不可に設定
    txt_b_31.config(state="normal")  # 一時的に編集可能に
    txt_b_31.delete(0, tk.END)
    txt_b_31.insert(0, round(UBtable['B'][2,0],4))
    txt_b_31.config(state="readonly") # 編集不可に設定
    txt_b_32.config(state="normal")  # 一時的に編集可能に
    txt_b_32.delete(0, tk.END)
    txt_b_32.insert(0, round(UBtable['B'][2,1],4))
    txt_b_32.config(state="readonly") # 編集不可に設定
    txt_b_33.config(state="normal")  # 一時的に編集可能に
    txt_b_33.delete(0, tk.END)
    txt_b_33.insert(0, round(UBtable['B'][2,2],4))
    txt_b_33.config(state="readonly") # 編集不可に設定
    
    # UB matrixを表示
    txt_ub_11.config(state="normal")  # 一時的に編集可能に
    txt_ub_11.delete(0, tk.END)
    txt_ub_11.insert(0, round(UBtable['UB'][0,0],4))
    txt_ub_11.config(state="readonly") # 編集不可に設定
    txt_ub_12.config(state="normal")  # 一時的に編集可能に
    txt_ub_12.delete(0, tk.END)
    txt_ub_12.insert(0, round(UBtable['UB'][0,1],4))
    txt_ub_12.config(state="readonly") # 編集不可に設定
    txt_ub_13.config(state="normal")  # 一時的に編集可能に
    txt_ub_13.delete(0, tk.END)
    txt_ub_13.insert(0, round(UBtable['UB'][0,2],4))
    txt_ub_13.config(state="readonly") # 編集不可に設定
    txt_ub_21.config(state="normal")  # 一時的に編集可能に
    txt_ub_21.delete(0, tk.END)
    txt_ub_21.insert(0, round(UBtable['UB'][1,0],4))
    txt_ub_21.config(state="readonly") # 編集不可に設定
    txt_ub_22.config(state="normal")  # 一時的に編集可能に
    txt_ub_22.delete(0, tk.END)
    txt_ub_22.insert(0, round(UBtable['UB'][1,1],4))
    txt_ub_22.config(state="readonly") # 編集不可に設定
    txt_ub_23.config(state="normal")  # 一時的に編集可能に
    txt_ub_23.delete(0, tk.END)
    txt_ub_23.insert(0, round(UBtable['UB'][1,2],4))
    txt_ub_23.config(state="readonly") # 編集不可に設定
    txt_ub_31.config(state="normal")  # 一時的に編集可能に
    txt_ub_31.delete(0, tk.END)
    txt_ub_31.insert(0, round(UBtable['UB'][2,0],4))
    txt_ub_31.config(state="readonly") # 編集不可に設定
    txt_ub_32.config(state="normal")  # 一時的に編集可能に
    txt_ub_32.delete(0, tk.END)
    txt_ub_32.insert(0, round(UBtable['UB'][2,1],4))
    txt_ub_32.config(state="readonly") # 編集不可に設定
    txt_ub_33.config(state="normal")  # 一時的に編集可能に
    txt_ub_33.delete(0, tk.END)
    txt_ub_33.insert(0, round(UBtable['UB'][2,2],4))
    txt_ub_33.config(state="readonly") # 編集不可に設定

# UB matrixの表示
"""
frame3s = ttk.Labelframe(frame3,text= "sense of spectrometer")
frame3s.grid(row=1,column=0,columnspan=2,sticky="NSEW")

frame3s.columnconfigure(0, weight=1)
frame3s.columnconfigure(1, weight=1)
frame3s.rowconfigure(0, weight=1)

# チェック有無変数
sense = tk.IntVar()
# value=0のラジオボタンにチェックを入れる
sense.set(0)

# ラジオボタンを設定
rdo_sense0 = tk.Radiobutton(frame3s, value=0, variable=sense, text='+++', width=15)
rdo_sense0.grid(row=0, column=0, sticky="NSEW")

rdo_sense1 = tk.Radiobutton(frame3s, value=1, variable=sense, text='-+-', width=15)
rdo_sense1.grid(row=0, column=1, sticky="NSEW")
"""

#ボタン1つで両方の計算を実行
UBcalculate_button = ttk.Button(frame3, text="display", command=displayUB)
UBcalculate_button.grid(row=1, column=1,sticky="NSEW")

# ブラッグピーク位置を入力
# ファイル選択のフレームの作成と設置
frame4 = ttk.Labelframe(root,text= "observed bragg peak position")
frame4.grid(row=3,column=0,sticky="NSEW")

frame4.columnconfigure(0, weight=1)
frame4.columnconfigure(1, weight=1)
frame4.columnconfigure(2, weight=1)
frame4.columnconfigure(3, weight=1)
frame4.columnconfigure(4, weight=1)
frame4.columnconfigure(5, weight=1)
frame4.columnconfigure(6, weight=1)
frame4.columnconfigure(7, weight=1)
frame4.rowconfigure(0, weight=1)
frame4.rowconfigure(1, weight=1)

def update_label():
    if eief.get() == 0:
        bpl1.config(text='Ei (meV)')
    else:
        bpl1.config(text='Ef (meV)')

# チェック有無変数
eief = tk.IntVar()
# value=0のラジオボタンにチェックを入れる
#eief.set(1)

# ラジオボタンを作成し、commandにupdate_labelを設定
rdo_eief0 = tk.Radiobutton(frame4, value=0, variable=eief, text='Ei fix', command=update_label, width=15)
rdo_eief0.grid(row=0, column=0, sticky="NSEW")

rdo_eief1 = tk.Radiobutton(frame4, value=1, variable=eief, text='Ef fix', command=update_label, width=15)
rdo_eief1.grid(row=1, column=0, sticky="NSEW")

if eief.get()==0:
    bpl1 = tk.Label(frame4,text='Ei (meV)')
elif eief.get()==1:
    bpl1 = tk.Label(frame4,text='Ef (meV)')
bpl1.grid(row=0, column=1,sticky="NSEW")

Energy = ttk.Entry(frame4)
Energy.grid(row=1, column=1,sticky="NSEW")

bpl2 = tk.Label(frame4,text='h')
bpl2.grid(row=0, column=2,sticky="NSEW")
bp_h = ttk.Entry(frame4)
bp_h.grid(row=1, column=2,sticky="NSEW")

bpl3 = tk.Label(frame4,text='k')
bpl3.grid(row=0, column=3,sticky="NSEW")
bp_k = ttk.Entry(frame4)
bp_k.grid(row=1, column=3,sticky="NSEW")

bpl4 = tk.Label(frame4,text='l')
bpl4.grid(row=0, column=4,sticky="NSEW")
bp_l = ttk.Entry(frame4)
bp_l.grid(row=1, column=4,sticky="NSEW")

bpl5 = tk.Label(frame4,text='C2')
bpl5.grid(row=0, column=5,sticky="NSEW")
bp_c2 = ttk.Entry(frame4)
bp_c2.grid(row=1, column=5,sticky="NSEW")

bpl6 = tk.Label(frame4,text='μ')
bpl6.grid(row=0, column=6,sticky="NSEW")
bp_mu = ttk.Entry(frame4)
bp_mu.grid(row=1, column=6,sticky="NSEW")

bpl7 = tk.Label(frame4,text='ν')
bpl7.grid(row=0, column=7,sticky="NSEW")
bp_nu = ttk.Entry(frame4)
bp_nu.grid(row=1, column=7,sticky="NSEW")

# select feature
# ファイル選択のフレームの作成と設置
frame6 = ttk.Labelframe(root,text= "select feature")
frame6.grid(row=4,column=0,sticky="NSEW")

# 計算する位置を入力
# テーマ
style = ttk.Style()

# タブの文字色を変えて見やすくする
# Style.map (TNotebook.Tab)
style.map(
    "example.TNotebook.Tab",
    foreground=[
        ('active', 'red'),
        ('disabled', 'black'),
        ('selected', 'blue'),
    ],
    background=[
        ('active', 'orange'),
        ('disabled', 'black'),
        ('selected', 'lightgreen'),
    ],
)

# Notebookウィジェットの作成
notebook00 = ttk.Notebook(frame6,style="example.TNotebook")

# タブの作成
tab_001 = tk.Frame(notebook00)# single Q-E point
tab_002 = tk.Frame(notebook00)# scan simulation
tab_003 = tk.Frame(notebook00)# lattice constant
tab_004 = tk.Frame(notebook00)# constant E scan
# notebookにタブを追加
notebook00.add(tab_001, text="single Q-E point")
notebook00.add(tab_002, text="scan simulation")
notebook00.add(tab_003, text="lattice constant")
notebook00.add(tab_004, text="toolbox")

# Notebookを配置
notebook00.pack(expand=True, fill="both")

# グリッドの重みを設定
tab_001.columnconfigure(0, weight=1)
tab_001.rowconfigure(0, weight=2)
tab_001.rowconfigure(1, weight=3)

tab_001a = ttk.Labelframe(tab_001,text= "calculation point")
tab_001a.grid(row=0,column=0,sticky="NSEW")
tab_001a.columnconfigure(0, weight=1)
tab_001a.columnconfigure(1, weight=1)
tab_001a.columnconfigure(2, weight=1)
tab_001a.columnconfigure(3, weight=1)
tab_001a.columnconfigure(4, weight=1)
tab_001a.rowconfigure(0, weight=1)
tab_001a.rowconfigure(1, weight=1)

def on_anglecalc():
    """angleを計算して返す"""
    params = get_parameters()
    
    # サンプル点の取得
    sv1 = np.array([float(sv1_h.get()), float(sv1_k.get()), float(sv1_l.get())])
    sv2 = np.array([float(sv2_h.get()), float(sv2_k.get()), float(sv2_l.get())])
    
    # RLtableを取得し、辞書から必要な変数を取り出す
    RLtable = on_Rlcalc()
    astar = RLtable['astar']
    bstar = RLtable['bstar']
    cstar = RLtable['cstar']
    alpha_star = RLtable['alpha_star']
    beta_star = RLtable['beta_star']
    gamma_star = RLtable['gamma_star']
    n_a = RLtable['n_a']
    n_b = RLtable['n_b']
    n_c = RLtable['n_c']
    
    # UBtableを取得し、辞書から必要な変数を取り出す
    UBtable = on_UBcalc()
    sv1 = UBtable['sv1']
    sv2 = UBtable['sv2']
    
    # ノルムの計算
    axis1 = sv1[0]*astar+sv1[1]*bstar+sv1[2]*cstar
    axis2 = sv2[0]*astar+sv2[1]*bstar+sv2[2]*cstar
    norm_sv1 = np.linalg.norm(sv1[0]*astar+sv1[1]*bstar+sv1[2]*cstar)
    norm_sv2 = np.linalg.norm(sv2[0]*astar+sv2[1]*bstar+sv2[2]*cstar)
    
    U=UBtable['U']
    B=UBtable['B']
    UB=UBtable['UB']
    
    # displayに表示
    # U matrixを表示
    txt_u_11.config(state="normal")  # 一時的に編集可能に
    txt_u_11.delete(0, tk.END)
    txt_u_11.insert(0, round(UBtable['U'][0,0],4))
    txt_u_11.config(state="readonly") # 編集不可に設定
    txt_u_12.config(state="normal")  # 一時的に編集可能に
    txt_u_12.delete(0, tk.END)
    txt_u_12.insert(0, round(UBtable['U'][0,1],4))
    txt_u_12.config(state="readonly") # 編集不可に設定
    txt_u_13.config(state="normal")  # 一時的に編集可能に
    txt_u_13.delete(0, tk.END)
    txt_u_13.insert(0, round(UBtable['U'][0,2],4))
    txt_u_13.config(state="readonly") # 編集不可に設定
    txt_u_21.config(state="normal")  # 一時的に編集可能に
    txt_u_21.delete(0, tk.END)
    txt_u_21.insert(0, round(UBtable['U'][1,0],4))
    txt_u_21.config(state="readonly") # 編集不可に設定
    txt_u_22.config(state="normal")  # 一時的に編集可能に
    txt_u_22.delete(0, tk.END)
    txt_u_22.insert(0, round(UBtable['U'][1,1],4))
    txt_u_22.config(state="readonly") # 編集不可に設定
    txt_u_23.config(state="normal")  # 一時的に編集可能に
    txt_u_23.delete(0, tk.END)
    txt_u_23.insert(0, round(UBtable['U'][1,2],4))
    txt_u_23.config(state="readonly") # 編集不可に設定
    txt_u_31.config(state="normal")  # 一時的に編集可能に
    txt_u_31.delete(0, tk.END)
    txt_u_31.insert(0, round(UBtable['U'][2,0],4))
    txt_u_31.config(state="readonly") # 編集不可に設定
    txt_u_32.config(state="normal")  # 一時的に編集可能に
    txt_u_32.delete(0, tk.END)
    txt_u_32.insert(0, round(UBtable['U'][2,1],4))
    txt_u_32.config(state="readonly") # 編集不可に設定
    txt_u_33.config(state="normal")  # 一時的に編集可能に
    txt_u_33.delete(0, tk.END)
    txt_u_33.insert(0, round(UBtable['U'][2,2],4))
    txt_u_33.config(state="readonly") # 編集不可に設定
    
    # B matrixを表示
    txt_b_11.config(state="normal")  # 一時的に編集可能に
    txt_b_11.delete(0, tk.END)
    txt_b_11.insert(0, round(UBtable['B'][0,0],4))
    txt_b_11.config(state="readonly") # 編集不可に設定
    txt_b_12.config(state="normal")  # 一時的に編集可能に
    txt_b_12.delete(0, tk.END)
    txt_b_12.insert(0, round(UBtable['B'][0,1],4))
    txt_b_12.config(state="readonly") # 編集不可に設定
    txt_b_13.config(state="normal")  # 一時的に編集可能に
    txt_b_13.delete(0, tk.END)
    txt_b_13.insert(0, round(UBtable['B'][0,2],4))
    txt_b_13.config(state="readonly") # 編集不可に設定
    txt_b_21.config(state="normal")  # 一時的に編集可能に
    txt_b_21.delete(0, tk.END)
    txt_b_21.insert(0, round(UBtable['B'][1,0],4))
    txt_b_21.config(state="readonly") # 編集不可に設定
    txt_b_22.config(state="normal")  # 一時的に編集可能に
    txt_b_22.delete(0, tk.END)
    txt_b_22.insert(0, round(UBtable['B'][1,1],4))
    txt_b_22.config(state="readonly") # 編集不可に設定
    txt_b_23.config(state="normal")  # 一時的に編集可能に
    txt_b_23.delete(0, tk.END)
    txt_b_23.insert(0, round(UBtable['B'][1,2],4))
    txt_b_23.config(state="readonly") # 編集不可に設定
    txt_b_31.config(state="normal")  # 一時的に編集可能に
    txt_b_31.delete(0, tk.END)
    txt_b_31.insert(0, round(UBtable['B'][2,0],4))
    txt_b_31.config(state="readonly") # 編集不可に設定
    txt_b_32.config(state="normal")  # 一時的に編集可能に
    txt_b_32.delete(0, tk.END)
    txt_b_32.insert(0, round(UBtable['B'][2,1],4))
    txt_b_32.config(state="readonly") # 編集不可に設定
    txt_b_33.config(state="normal")  # 一時的に編集可能に
    txt_b_33.delete(0, tk.END)
    txt_b_33.insert(0, round(UBtable['B'][2,2],4))
    txt_b_33.config(state="readonly") # 編集不可に設定
    
    # UB matrixを表示
    txt_ub_11.config(state="normal")  # 一時的に編集可能に
    txt_ub_11.delete(0, tk.END)
    txt_ub_11.insert(0, round(UBtable['UB'][0,0],4))
    txt_ub_11.config(state="readonly") # 編集不可に設定
    txt_ub_12.config(state="normal")  # 一時的に編集可能に
    txt_ub_12.delete(0, tk.END)
    txt_ub_12.insert(0, round(UBtable['UB'][0,1],4))
    txt_ub_12.config(state="readonly") # 編集不可に設定
    txt_ub_13.config(state="normal")  # 一時的に編集可能に
    txt_ub_13.delete(0, tk.END)
    txt_ub_13.insert(0, round(UBtable['UB'][0,2],4))
    txt_ub_13.config(state="readonly") # 編集不可に設定
    txt_ub_21.config(state="normal")  # 一時的に編集可能に
    txt_ub_21.delete(0, tk.END)
    txt_ub_21.insert(0, round(UBtable['UB'][1,0],4))
    txt_ub_21.config(state="readonly") # 編集不可に設定
    txt_ub_22.config(state="normal")  # 一時的に編集可能に
    txt_ub_22.delete(0, tk.END)
    txt_ub_22.insert(0, round(UBtable['UB'][1,1],4))
    txt_ub_22.config(state="readonly") # 編集不可に設定
    txt_ub_23.config(state="normal")  # 一時的に編集可能に
    txt_ub_23.delete(0, tk.END)
    txt_ub_23.insert(0, round(UBtable['UB'][1,2],4))
    txt_ub_23.config(state="readonly") # 編集不可に設定
    txt_ub_31.config(state="normal")  # 一時的に編集可能に
    txt_ub_31.delete(0, tk.END)
    txt_ub_31.insert(0, round(UBtable['UB'][2,0],4))
    txt_ub_31.config(state="readonly") # 編集不可に設定
    txt_ub_32.config(state="normal")  # 一時的に編集可能に
    txt_ub_32.delete(0, tk.END)
    txt_ub_32.insert(0, round(UBtable['UB'][2,1],4))
    txt_ub_32.config(state="readonly") # 編集不可に設定
    txt_ub_33.config(state="normal")  # 一時的に編集可能に
    txt_ub_33.delete(0, tk.END)
    txt_ub_33.insert(0, round(UBtable['UB'][2,2],4))
    txt_ub_33.config(state="readonly") # 編集不可に設定
    
    # Bragg peak positionの取得
    bpe = float(Energy.get())
    bpc2 = float(bp_c2.get())
    bpmu = float(bp_mu.get())
    bpnu = float(bp_nu.get())
    #bp = np.array([float(bp_h.get()), float(bp_k.get()), float(bp_l.get())])
    
    bp_h_entry = bp_h.get()
    try:
        # 少数の場合
        if '/' not in bp_h_entry:
            bph = float(bp_h_entry)
        else:
            # 分数の場合
            bph = float(Fraction(bp_h_entry))
    except ValueError:
        pass
    
    bp_k_entry = bp_k.get()
    try:
        # 少数の場合
        if '/' not in bp_k_entry:
            bpk = float(bp_k_entry)
        else:
            # 分数の場合
            bpk = float(Fraction(bp_k_entry))
    except ValueError:
        pass
    
    bp_l_entry = bp_l.get()
    try:
        # 少数の場合
        if '/' not in bp_l_entry:
            bpl = float(bp_l_entry)
        else:
            # 分数の場合
            bpl = float(Fraction(bp_l_entry))
    except ValueError:
        pass
    
    bp = np.array([bph,bpk,bpl])
    
    # calculation pointの取得
    cphw = float(acbe.get())
    #cp = np.array([float(acbh.get()), float(acbk.get()), float(acbl.get())])
    
    cp_h_entry = acbh.get()
    try:
        # 少数の場合
        if '/' not in cp_h_entry:
            cph = float(cp_h_entry)
        else:
            # 分数の場合
            cph = float(Fraction(cp_h_entry))
    except ValueError:
        pass
    
    cp_k_entry = acbk.get()
    try:
        # 少数の場合
        if '/' not in cp_k_entry:
            cpk = float(cp_k_entry)
        else:
            # 分数の場合
            cpk = float(Fraction(cp_k_entry))
    except ValueError:
        pass
    
    cp_l_entry = acbl.get()
    try:
        # 少数の場合
        if '/' not in cp_l_entry:
            cpl = float(cp_l_entry)
        else:
            # 分数の場合
            cpl = float(Fraction(cp_l_entry))
    except ValueError:
        pass
    
    cp = np.array([cph,cpk,cpl])
    
    # Ei or Ef fixの判定
    fixe=float(eief.get())
    
    # angleを計算
    angletable = angle_calc(
        astar,bstar,cstar,U,B,UB,bpe,bpc2,bpmu,bpnu,bp,cphw,cp,fixe
    )
    
    return angletable

def calculate_angle():
    # on_anglecalc の結果を取得
    angletable,error_message = on_anglecalc()
    #print("RLtable:", RLtable)
    
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
    
    monochromator_radius = float(config['settings']['monochromator_radius'])
    monochromator_to_sample = float(config['settings']['monochromator_to_sample'])
    sample_goniometer_radius = float(config['settings']['sample_goniometer_radius'])
    sample_to_analyzer = float(config['settings']['sample_to_analyzer'])
    analyzer_radius = float(config['settings']['analyzer_radius'])
    analyzer_to_detector = float(config['settings']['analyzer_to_detector'])
    detector_radius = float(config['settings']['detector_radius'])
    floor_length = float(config['settings']['floor_length'])
    floor_width = float(config['settings']['floor_width'])
    floor_position_x = float(config['settings']['floor_position_x'])
    floor_position_y = float(config['settings']['floor_position_y'])
    
    # sample gonioがfloorからはみ出る場合
    positionY_sample = monochromator_to_sample * np.sin(np.radians(angletable['A1'])) - sample_goniometer_radius # < floor_position_y
    # analyzer dramがfloorからはみ出る場合
    positionY_analyzer = monochromator_to_sample * np.sin(np.radians(angletable['A1'])) - sample_to_analyzer * np.sin(np.radians(angletable['A2'] - angletable['A1'])) - analyzer_radius # < floor_position_y
    # detector dramがfloorからはみ出る場合
    positionY_detector = monochromator_to_sample * np.sin(np.radians(angletable['A1'])) - sample_to_analyzer * np.sin(np.radians(angletable['A2'] - angletable['A1'])) + analyzer_to_detector * np.sin(np.radians(angletable['A3'] - angletable['A2'] + angletable['A1'])) - detector_radius # < floor_position_y
    positionX_detector = monochromator_to_sample * np.cos(np.radians(angletable['A1'])) + sample_to_analyzer * np.cos(np.radians(angletable['A2'] - angletable['A1'])) + analyzer_to_detector * np.cos(np.radians(angletable['A3'] - angletable['A2'] + angletable['A1'])) - detector_radius # < floor_position_x+floor_length
    
    if error_message is not None:
        acl10.config(text=error_message, fg="red")
        return  # 計算を中断
    elif (positionY_sample < floor_position_y or
          positionY_analyzer < floor_position_y or
          positionY_detector < floor_position_y or
          positionX_detector > floor_position_x+floor_length):
        acl10.config(text="Out of the floor.", fg="red")
    elif (round(angletable['C1'],4)<float(hwl2f.get()) or
        round(angletable['C1'],4)>float(hwl2t.get()) or
        round(angletable['A1'],4)<float(hwl3f.get()) or
        round(angletable['A1'],4)>float(hwl3t.get()) or
        round(angletable['C2'],4)<float(hwl4f.get()) or
        round(angletable['C2'],4)>float(hwl4t.get()) or
        round(angletable['A2'],4)<float(hwl5f.get()) or
        round(angletable['A2'],4)>float(hwl5t.get()) or
        round(angletable['C3'],4)<float(hwl6f.get()) or
        round(angletable['C3'],4)>float(hwl6t.get()) or
        round(angletable['A3'],4)<float(hwl7f.get()) or
        round(angletable['A3'],4)>float(hwl7t.get()) or
        round(angletable['mu'],4)<float(hwl8f.get()) or
        round(angletable['mu'],4)>float(hwl8t.get()) or
        round(angletable['nu'],4)<float(hwl9f.get()) or
        round(angletable['nu'],4)>float(hwl9t.get())):
        acl10.config(text="Out of the hardware limit range.", fg="red")
    
    else:
        # `phi_cal` の結果を用いた続きの処理
        acl10.config(text="The calculation was successful.", fg="black")
    
    # angle計算を表示
    acb1.config(state="normal")  # 一時的に編集可能に
    acb1.delete(0, tk.END)
    acb1.insert(0, round(angletable['C1'],4))
    acb1.config(state="readonly") # 編集不可に設定
    
    acb2.config(state="normal")  # 一時的に編集可能に
    acb2.delete(0, tk.END)
    acb2.insert(0, round(angletable['A1'],4))
    acb2.config(state="readonly") # 編集不可に設定
    
    acb3.config(state="normal")  # 一時的に編集可能に
    acb3.delete(0, tk.END)
    acb3.insert(0, round(angletable['C2'],4))
    acb3.config(state="readonly") # 編集不可に設定
    
    acb4.config(state="normal")  # 一時的に編集可能に
    acb4.delete(0, tk.END)
    acb4.insert(0, round(angletable['A2'],4))
    acb4.config(state="readonly") # 編集不可に設定
    
    acb5.config(state="normal")  # 一時的に編集可能に
    acb5.delete(0, tk.END)
    acb5.insert(0, round(angletable['C3'],4))
    acb5.config(state="readonly") # 編集不可に設定
    
    acb6.config(state="normal")  # 一時的に編集可能に
    acb6.delete(0, tk.END)
    acb6.insert(0, round(angletable['A3'],4))
    acb6.config(state="readonly") # 編集不可に設定
    
    acb7.config(state="normal")  # 一時的に編集可能に
    acb7.delete(0, tk.END)
    acb7.insert(0, round(angletable['mu'],4))
    acb7.config(state="readonly") # 編集不可に設定
    
    acb8.config(state="normal")  # 一時的に編集可能に
    acb8.delete(0, tk.END)
    acb8.insert(0, round(angletable['nu'],4))
    acb8.config(state="readonly") # 編集不可に設定
    
    A_sets = []  # A_setsリストを初期化
    QE_sets = []
    
    # A1, A2, A3 を取得して A_sets に追加
    A1 = round(angletable['A1'],4)  # 'A1'
    A2 = -round(angletable['A2'],4)  # 'A2'
    A3 = round(angletable['A3'],4)  # 'A3'
    A_sets.append([A1, A2, A3])  # A_sets に追加
    # hw, h,k,l
    hw = round(float(acbe.get()), 4)  # 'hw'
    #h = round(float(acbh.get()), 4)  # 'h'
    #k = round(float(acbk.get()), 4)  # 'k'
    #l = round(float(acbl.get()), 4)  # 'l'
    
    cp_h_entry = acbh.get()
    try:
        # 少数の場合
        if '/' not in cp_h_entry:
            h = float(cp_h_entry)
        else:
            # 分数の場合
            h = float(Fraction(cp_h_entry))
    except ValueError:
        pass
    
    cp_k_entry = acbk.get()
    try:
        # 少数の場合
        if '/' not in cp_k_entry:
            k = float(cp_k_entry)
        else:
            # 分数の場合
            k = float(Fraction(cp_k_entry))
    except ValueError:
        pass
    
    cp_l_entry = acbl.get()
    try:
        # 少数の場合
        if '/' not in cp_l_entry:
            l = float(cp_l_entry)
        else:
            # 分数の場合
            l = float(Fraction(cp_l_entry))
    except ValueError:
        pass
    
    h=round(h,4)
    k=round(k,4)
    l=round(l,4)
    
    QE_sets.append([hw, h, k, l])

    if fig_spec.get()==1:
        # プロット関数を呼び出し
        #plot_spectrometer_with_gif(A_sets,QE_sets)
        plot_spectrometer(A_sets,QE_sets)
    
    Hfocus = calc_hf.get()
    num_ana = float(acna.get())
    
    if fig_reci.get()==1:
        # 逆格子空間のki,kf,τベクトルを示す。
        params = get_parameters()
        RLtable = RL_calc(**params)  # RL_calcに辞書を展開して渡す
        
        # 散乱面の取得
        sv1 = np.array([float(sv1_h.get()), float(sv1_k.get()), float(sv1_l.get())])
        sv2 = np.array([float(sv2_h.get()), float(sv2_k.get()), float(sv2_l.get())])
        
        cp_h_entry = acbh.get()
        try:
            # 少数の場合
            if '/' not in cp_h_entry:
                cph = float(cp_h_entry)
            else:
                # 分数の場合
                cph = float(Fraction(cp_h_entry))
        except ValueError:
            pass
        
        cp_k_entry = acbk.get()
        try:
            # 少数の場合
            if '/' not in cp_k_entry:
                cpk = float(cp_k_entry)
            else:
                # 分数の場合
                cpk = float(Fraction(cp_k_entry))
        except ValueError:
            pass
        
        cp_l_entry = acbl.get()
        try:
            # 少数の場合
            if '/' not in cp_l_entry:
                cpl = float(cp_l_entry)
            else:
                # 分数の場合
                cpl = float(Fraction(cp_l_entry))
        except ValueError:
            pass
        
        cphw = float(acbe.get())
        cp = np.array([cph,cpk,cpl])
        
        fixe=float(eief.get())
        bpe = float(Energy.get())
        
        plot_reciprocal_space(bpe,cphw,cp,fixe,sv1,sv2,**RLtable)
    
    if fig_reso.get()==1:
        bpe = float(Energy.get())
        fixe=float(eief.get())
        
        # Entry ウィジェットの値を辞書にまとめる
        entry_values = {
            "div_1st_m": div_1st_m.get(),
            "div_2nd_h": div_2nd_h.get(),
            "div_2nd_v": div_2nd_v.get(),
            "div_3rd_h": div_3rd_h.get(),
            "div_3rd_v": div_3rd_v.get(),
            "div_4th_h": div_4th_h.get(),
            "div_4th_v": div_4th_v.get(),
            "mos_mono_h": mos_mono_h.get(),
            "mos_mono_v": mos_mono_v.get(),
            "mos_sam_h": mos_sam_h.get(),
            "mos_sam_v": mos_sam_v.get(),
            "mos_ana_h": mos_ana_h.get(),
            "mos_ana_v": mos_ana_v.get(),
        }
        calcresolution(A_sets,QE_sets,bpe,fixe,hw,Hfocus,num_ana,entry_values)
    
    plt.show()

acle = tk.Label(tab_001a,text='ℏω')
acle.grid(row=0, column=0,sticky="NSEW")
acbe = ttk.Entry(tab_001a)
acbe.grid(row=1, column=0,sticky="NSEW")
acbe.insert(0,'0')

aclh = tk.Label(tab_001a,text='h')
aclh.grid(row=0, column=1,sticky="NSEW")
acbh = ttk.Entry(tab_001a)
acbh.grid(row=1, column=1,sticky="NSEW")
acbh.insert(0,'1')

aclk = tk.Label(tab_001a,text='k')
aclk.grid(row=0, column=2,sticky="NSEW")
acbk = ttk.Entry(tab_001a)
acbk.grid(row=1, column=2,sticky="NSEW")
acbk.insert(0,'0')

acll = tk.Label(tab_001a,text='l')
acll.grid(row=0, column=3,sticky="NSEW")
acbl = ttk.Entry(tab_001a)
acbl.grid(row=1, column=3,sticky="NSEW")
acbl.insert(0,'0')

#ボタン1つで両方の計算を実行
Angle_calculate_button = ttk.Button(tab_001a, text="calc", command=calculate_angle,width=16)
Angle_calculate_button.grid(row=1, column=4,sticky="NSEW")

tab_001b = ttk.Labelframe(tab_001,text= "calculation results (unit : deg)")
tab_001b.grid(row=1,column=0,sticky="NSEW")
tab_001b.columnconfigure(0, weight=1)
tab_001b.columnconfigure(1, weight=1)
tab_001b.columnconfigure(2, weight=1)
tab_001b.columnconfigure(3, weight=1)
tab_001b.columnconfigure(4, weight=1)
tab_001b.columnconfigure(5, weight=1)
tab_001b.columnconfigure(6, weight=1)
tab_001b.columnconfigure(7, weight=1)
tab_001b.rowconfigure(0, weight=1)
tab_001b.rowconfigure(1, weight=1)
tab_001b.rowconfigure(2, weight=1)

acl1 = tk.Label(tab_001b,text='C1')
acl1.grid(row=0, column=0,sticky="NSEW")
acb1 = ttk.Entry(tab_001b)
acb1.grid(row=1, column=0,sticky="NSEW")
acb1.insert(0,'0')
acb1.config(state="readonly") # 編集不可に設定

acl2 = tk.Label(tab_001b,text='A1')
acl2.grid(row=0, column=1,sticky="NSEW")
acb2 = ttk.Entry(tab_001b)
acb2.grid(row=1, column=1,sticky="NSEW")
acb2.insert(0,'0')
acb2.config(state="readonly") # 編集不可に設定

acl3 = tk.Label(tab_001b,text='C2')
acl3.grid(row=0, column=2,sticky="NSEW")
acb3 = ttk.Entry(tab_001b)
acb3.grid(row=1, column=2,sticky="NSEW")
acb3.insert(0,'0')
acb3.config(state="readonly") # 編集不可に設定

acl4 = tk.Label(tab_001b,text='A2')
acl4.grid(row=0, column=3,sticky="NSEW")
acb4 = ttk.Entry(tab_001b)
acb4.grid(row=1, column=3,sticky="NSEW")
acb4.insert(0,'0')
acb4.config(state="readonly") # 編集不可に設定

acl5 = tk.Label(tab_001b,text='C3')
acl5.grid(row=0, column=4,sticky="NSEW")
acb5 = ttk.Entry(tab_001b)
acb5.grid(row=1, column=4,sticky="NSEW")
acb5.insert(0,'0')
acb5.config(state="readonly") # 編集不可に設定

acl6 = tk.Label(tab_001b,text='A3')
acl6.grid(row=0, column=5,sticky="NSEW")
acb6 = ttk.Entry(tab_001b)
acb6.grid(row=1, column=5,sticky="NSEW")
acb6.insert(0,'0')
acb6.config(state="readonly") # 編集不可に設定

acl7 = tk.Label(tab_001b,text='μ')
acl7.grid(row=0, column=6,sticky="NSEW")
acb7 = ttk.Entry(tab_001b)
acb7.grid(row=1, column=6,sticky="NSEW")
acb7.insert(0,'0')
acb7.config(state="readonly") # 編集不可に設定

acl8 = tk.Label(tab_001b,text='ν')
acl8.grid(row=0, column=7,sticky="NSEW")
acb8 = ttk.Entry(tab_001b)
acb8.grid(row=1, column=7,sticky="NSEW")
acb8.insert(0,'0')
acb8.config(state="readonly") # 編集不可に設定

acl9 = tk.Label(tab_001b,text='warning : ')
acl9.grid(row=2, column=0,sticky="NSEW")

acl10 = tk.Label(tab_001b,text='')
acl10.grid(row=2, column=1,columnspan=7,sticky="NSEW")

# collimator information
# ファイル選択のフレームの作成と設置
frame6 = ttk.Labelframe(root,text= "instrumental condition")
frame6.grid(row=5,column=0,sticky="NSEW")

frame6.columnconfigure(0, weight=5)
frame6.columnconfigure(1, weight=4)
frame6.columnconfigure(2, weight=1)
frame6.rowconfigure(0, weight=1)

frame6a = ttk.Labelframe(frame6,text= "divergence (unit : min)")
frame6a.grid(row=0,column=0,sticky="NSEW")

frame6a.columnconfigure(0, weight=1)
frame6a.columnconfigure(1, weight=1)
frame6a.columnconfigure(2, weight=1)
frame6a.columnconfigure(3, weight=1)
frame6a.columnconfigure(4, weight=1)
frame6a.rowconfigure(0, weight=1)
frame6a.rowconfigure(1, weight=1)
frame6a.rowconfigure(2, weight=1)

label1 = tk.Label(frame6a,text='H',width = 12)
label1.grid(row=1, column=1,sticky="NSEW")
label2 = tk.Label(frame6a,text='V',width = 12)
label2.grid(row=2, column=1,sticky="NSEW")
label3 = tk.Label(frame6a,text='guide')
label3.grid(row=0, column=0,sticky="NSEW")
label3m = tk.Label(frame6a,text='m')
label3m.grid(row=1, column=0,sticky="NSEW")
label4 = tk.Label(frame6a,text='2nd col')
label4.grid(row=0, column=2,sticky="NSEW")
label5 = tk.Label(frame6a,text='3rd col')
label5.grid(row=0, column=3,sticky="NSEW")
label6 = tk.Label(frame6a,text='4th col')
label6.grid(row=0, column=4,sticky="NSEW")

div_1st_m = ttk.Entry(frame6a)
div_1st_m.grid(row=2, column=0,sticky="NSEW")
div_2nd_h = ttk.Entry(frame6a)
div_2nd_h.grid(row=1, column=2,sticky="NSEW")
div_2nd_v = ttk.Entry(frame6a)
div_2nd_v.grid(row=2, column=2,sticky="NSEW")
div_3rd_h = ttk.Entry(frame6a)
div_3rd_h.grid(row=1, column=3,sticky="NSEW")
div_3rd_v = ttk.Entry(frame6a)
div_3rd_v.grid(row=2, column=3,sticky="NSEW")
div_4th_h = ttk.Entry(frame6a)
div_4th_h.grid(row=1, column=4,sticky="NSEW")
div_4th_v = ttk.Entry(frame6a)
div_4th_v.grid(row=2, column=4,sticky="NSEW")

frame6b = ttk.Labelframe(frame6,text= "mosaic (unit : min)")
frame6b.grid(row=0,column=1,sticky="NSEW")

frame6b.columnconfigure(0, weight=1)
frame6b.columnconfigure(1, weight=1)
frame6b.columnconfigure(2, weight=1)
frame6b.columnconfigure(3, weight=1)
frame6b.rowconfigure(0, weight=1)
frame6b.rowconfigure(1, weight=1)
frame6b.rowconfigure(2, weight=1)

label1 = tk.Label(frame6b,text='H',width = 12)
label1.grid(row=1, column=0,sticky="NSEW")
label2 = tk.Label(frame6b,text='V',width = 12)
label2.grid(row=2, column=0,sticky="NSEW")
label7 = tk.Label(frame6b,text='monochro')
label7.grid(row=0, column=1,sticky="NSEW")
label8 = tk.Label(frame6b,text='sample')
label8.grid(row=0, column=2,sticky="NSEW")
label8 = tk.Label(frame6b,text='analyzer')
label8.grid(row=0, column=3,sticky="NSEW")

mos_mono_h = ttk.Entry(frame6b)
mos_mono_h.grid(row=1, column=1,sticky="NSEW")
mos_mono_v = ttk.Entry(frame6b)
mos_mono_v.grid(row=2, column=1,sticky="NSEW")

mos_sam_h = ttk.Entry(frame6b)
mos_sam_h.grid(row=1, column=2,sticky="NSEW")
mos_sam_v = ttk.Entry(frame6b)
mos_sam_v.grid(row=2, column=2,sticky="NSEW")

mos_ana_h = ttk.Entry(frame6b)
mos_ana_h.grid(row=1, column=3,sticky="NSEW")
mos_ana_v = ttk.Entry(frame6b)
mos_ana_v.grid(row=2, column=3,sticky="NSEW")

frame6c = ttk.Labelframe(frame6,text= "focusing")
frame6c.grid(row=0,column=2,sticky="NSEW")
frame6c.columnconfigure(0, weight=1)
frame6c.rowconfigure(0, weight=1)
frame6c.rowconfigure(1, weight=1)
frame6c.rowconfigure(2, weight=1)

#チェック有無変数
calc_hf = tk.IntVar()
# value=0にチェックを入れる
calc_hf.set(0)

calc_HF = tk.Checkbutton(frame6c, variable=calc_hf, text='on')
calc_HF.grid(row=0, column=0,sticky="NSEW")

aclna = tk.Label(frame6c,text='blade num')
aclna.grid(row=1, column=0,sticky="NSEW")
acna = ttk.Entry(frame6c)
acna.grid(row=2, column=0,sticky="NSEW")
acna.insert(0,'7')

# ファイル選択のフレームの作成と設置
frame5 = ttk.Labelframe(root)
frame5.grid(row=6,column=0,sticky="NSEW")

frame5.columnconfigure(0, weight=9)
frame5.columnconfigure(1, weight=1)
frame5.rowconfigure(0, weight=1)

# hardware limit
# ファイル選択のフレームの作成と設置
frame5a = ttk.Labelframe(frame5,text= "hardware limit")
frame5a.grid(row=0,column=0,sticky="NSEW")

frame5a.columnconfigure(0, weight=1)
frame5a.columnconfigure(1, weight=1)
frame5a.columnconfigure(2, weight=1)
frame5a.columnconfigure(3, weight=1)
frame5a.columnconfigure(4, weight=1)
frame5a.columnconfigure(5, weight=1)
frame5a.columnconfigure(6, weight=1)
frame5a.columnconfigure(7, weight=1)
frame5a.columnconfigure(8, weight=1)
frame5a.rowconfigure(0, weight=1)
frame5a.rowconfigure(1, weight=1)
frame5a.rowconfigure(2, weight=1)

hwll0 = tk.Label(frame5a,text='from',width=15)
hwll0.grid(row=1, column=0,sticky="NSEW")
hwll1 = tk.Label(frame5a,text='to',width=15)
hwll1.grid(row=2, column=0,sticky="NSEW")

hwll2 = tk.Label(frame5a,text='C1')
hwll2.grid(row=0, column=1,sticky="NSEW")
hwl2f = ttk.Entry(frame5a)
hwl2f.grid(row=1, column=1,sticky="NSEW")
#hwl2f.insert(0,'19.9305')
hwl2t = ttk.Entry(frame5a)
hwl2t.grid(row=2, column=1,sticky="NSEW")
#hwl2t.insert(0,'58.482')

hwll3 = tk.Label(frame5a,text='A1')
hwll3.grid(row=0, column=2,sticky="NSEW")
hwl3f = ttk.Entry(frame5a)
hwl3f.grid(row=1, column=2,sticky="NSEW")
#hwl3f.insert(0,'39.861')
hwl3t = ttk.Entry(frame5a)
hwl3t.grid(row=2, column=2,sticky="NSEW")
#hwl3t.insert(0,'116.964')

hwll4 = tk.Label(frame5a,text='C2')
hwll4.grid(row=0, column=3,sticky="NSEW")
hwl4f = ttk.Entry(frame5a)
hwl4f.grid(row=1, column=3,sticky="NSEW")
#hwl4f.insert(0,'-180')
hwl4t = ttk.Entry(frame5a)
hwl4t.grid(row=2, column=3,sticky="NSEW")
#hwl4t.insert(0,'180')

hwll5 = tk.Label(frame5a,text='A2')
hwll5.grid(row=0, column=4,sticky="NSEW")
hwl5f = ttk.Entry(frame5a)
hwl5f.grid(row=1, column=4,sticky="NSEW")
#hwl5f.insert(0,'6')
hwl5t = ttk.Entry(frame5a)
hwl5t.grid(row=2, column=4,sticky="NSEW")
#hwl5t.insert(0,'120')

hwll6 = tk.Label(frame5a,text='C3')
hwll6.grid(row=0, column=5,sticky="NSEW")
hwl6f = ttk.Entry(frame5a)
hwl6f.grid(row=1, column=5,sticky="NSEW")
#hwl6f.insert(0,'19.9305')
hwl6t = ttk.Entry(frame5a)
hwl6t.grid(row=2, column=5,sticky="NSEW")
#hwl6t.insert(0,'58.482')

hwll7 = tk.Label(frame5a,text='A3')
hwll7.grid(row=0, column=6,sticky="NSEW")
hwl7f = ttk.Entry(frame5a)
hwl7f.grid(row=1, column=6,sticky="NSEW")
#hwl7f.insert(0,'39.861')
hwl7t = ttk.Entry(frame5a)
hwl7t.grid(row=2, column=6,sticky="NSEW")
#hwl7t.insert(0,'116.964')

hwll8 = tk.Label(frame5a,text='μ')
hwll8.grid(row=0, column=7,sticky="NSEW")
hwl8f = ttk.Entry(frame5a)
hwl8f.grid(row=1, column=7,sticky="NSEW")
#hwl8f.insert(0,'-5')
hwl8t = ttk.Entry(frame5a)
hwl8t.grid(row=2, column=7,sticky="NSEW")
#hwl8t.insert(0,'5')

hwll9 = tk.Label(frame5a,text='ν')
hwll9.grid(row=0, column=8,sticky="NSEW")
hwl9f = ttk.Entry(frame5a)
hwl9f.grid(row=1, column=8,sticky="NSEW")
#hwl9f.insert(0,'-5')
hwl9t = ttk.Entry(frame5a)
hwl9t.grid(row=2, column=8,sticky="NSEW")
#hwl9t.insert(0,'5')

frame5b = ttk.Labelframe(frame5,text= "figure")
frame5b.grid(row=0,column=1,sticky="NSEW")

frame5b.columnconfigure(0, weight=1)
frame5b.rowconfigure(0, weight=1)
frame5b.rowconfigure(1, weight=1)
frame5b.rowconfigure(2, weight=1)

# 分光器図と分解能を表示するかどうかのチェックボックス
# チェック有無変数
fig_spec = tk.IntVar()
# value=0にチェックを入れる
fig_spec.set(1)

show_fig_spec = tk.Checkbutton(frame5b, variable=fig_spec, text='spec',width=14)
show_fig_spec.grid(row=0, column=0,sticky="NSEW")

# チェック有無変数
fig_reci = tk.IntVar()
# value=0にチェックを入れる
fig_reci.set(1)

show_fig_reci = tk.Checkbutton(frame5b, variable=fig_reci, text='reci')
show_fig_reci.grid(row=1, column=0,sticky="NSEW")

# チェック有無変数
fig_reso = tk.IntVar()
# value=0にチェックを入れる
fig_reso.set(1)

show_fig_reso = tk.Checkbutton(frame5b, variable=fig_reso, text='reso')
show_fig_reso.grid(row=2, column=0,sticky="NSEW")

# グリッドの重みを設定
tab_002.columnconfigure(0, weight=1)
tab_002.columnconfigure(1, weight=1)
tab_002.rowconfigure(0, weight=1)

tab_002a = ttk.Labelframe(tab_002,text= "constant Q scan")
tab_002a.grid(row=0,column=0,sticky="NSEW")
tab_002a.columnconfigure(0, weight=1)
tab_002a.columnconfigure(1, weight=1)
tab_002a.columnconfigure(2, weight=1)
tab_002a.columnconfigure(3, weight=1)
tab_002a.columnconfigure(4, weight=1)
tab_002a.rowconfigure(0, weight=1)
tab_002a.rowconfigure(1, weight=1)
tab_002a.rowconfigure(2, weight=1)
tab_002a.rowconfigure(3, weight=1)
tab_002a.rowconfigure(4, weight=1)

cqslt = tk.Label(tab_002a,text='ℏω')
cqslt.grid(row=0, column=1,sticky="NSEW")
cqslf = tk.Label(tab_002a,text='from',width=15)
cqslf.grid(row=1, column=0,sticky="NSEW")
cqslt = tk.Label(tab_002a,text='to')
cqslt.grid(row=2, column=0,sticky="NSEW")
cqsli = tk.Label(tab_002a,text='inc')
cqsli.grid(row=3, column=0,sticky="NSEW")

cqsl1 = tk.Label(tab_002a,text='h')
cqsl1.grid(row=0, column=2,sticky="NSEW")
cqsl2 = tk.Label(tab_002a,text='k')
cqsl2.grid(row=0, column=3,sticky="NSEW")
cqsl2 = tk.Label(tab_002a,text='l')
cqsl2.grid(row=0, column=4,sticky="NSEW")

cqsef = ttk.Entry(tab_002a)
cqsef.grid(row=1, column=1,sticky="NSEW")
cqsef.insert(0,'-0.4')
cqset = ttk.Entry(tab_002a)
cqset.grid(row=2, column=1,sticky="NSEW")
cqset.insert(0,'7')
cqsei = ttk.Entry(tab_002a)
cqsei.grid(row=3, column=1,sticky="NSEW")
cqsei.insert(0,'0.2')

cqse1 = ttk.Entry(tab_002a)
cqse1.grid(row=1, column=2,sticky="NSEW")
cqse1.insert(0,'1.5')
cqse2 = ttk.Entry(tab_002a)
cqse2.grid(row=1, column=3,sticky="NSEW")
cqse2.insert(0,'0')
cqse3 = ttk.Entry(tab_002a)
cqse3.grid(row=1, column=4,sticky="NSEW")
cqse3.insert(0,'0')

# 別ウィンドウでデータテーブルを表示する関数
def constQscan_show_table():
    # 新しいウィンドウを作成
    result_window = tk.Toplevel()
    result_window.title("calculation results (unit : deg)")
    
    # Treeviewの設定
    tree = ttk.Treeview(result_window, columns=("hw","h","k","l","C1", "A1", "C2", "A2", "C3", "A3", "mu", "nu"), show="headings")
    tree.pack(fill="both", expand=True)
    
    # 各列に見出しを設定
    for col in tree["columns"]:
        tree.heading(col, text=col)
        tree.column(col, width=80, anchor="center")
    
    # RLtableを取得し、辞書から必要な変数を取り出す
    RLtable = on_Rlcalc()
    astar = RLtable['astar']
    bstar = RLtable['bstar']
    cstar = RLtable['cstar']
    
    # UBtableを取得し、辞書から必要な変数を取り出す
    UBtable = on_UBcalc()
    
    U=UBtable['U']
    B=UBtable['B']
    UB=UBtable['UB']
    
    # displayに表示
    # U matrixを表示
    txt_u_11.config(state="normal")  # 一時的に編集可能に
    txt_u_11.delete(0, tk.END)
    txt_u_11.insert(0, round(UBtable['U'][0,0],4))
    txt_u_11.config(state="readonly") # 編集不可に設定
    txt_u_12.config(state="normal")  # 一時的に編集可能に
    txt_u_12.delete(0, tk.END)
    txt_u_12.insert(0, round(UBtable['U'][0,1],4))
    txt_u_12.config(state="readonly") # 編集不可に設定
    txt_u_13.config(state="normal")  # 一時的に編集可能に
    txt_u_13.delete(0, tk.END)
    txt_u_13.insert(0, round(UBtable['U'][0,2],4))
    txt_u_13.config(state="readonly") # 編集不可に設定
    txt_u_21.config(state="normal")  # 一時的に編集可能に
    txt_u_21.delete(0, tk.END)
    txt_u_21.insert(0, round(UBtable['U'][1,0],4))
    txt_u_21.config(state="readonly") # 編集不可に設定
    txt_u_22.config(state="normal")  # 一時的に編集可能に
    txt_u_22.delete(0, tk.END)
    txt_u_22.insert(0, round(UBtable['U'][1,1],4))
    txt_u_22.config(state="readonly") # 編集不可に設定
    txt_u_23.config(state="normal")  # 一時的に編集可能に
    txt_u_23.delete(0, tk.END)
    txt_u_23.insert(0, round(UBtable['U'][1,2],4))
    txt_u_23.config(state="readonly") # 編集不可に設定
    txt_u_31.config(state="normal")  # 一時的に編集可能に
    txt_u_31.delete(0, tk.END)
    txt_u_31.insert(0, round(UBtable['U'][2,0],4))
    txt_u_31.config(state="readonly") # 編集不可に設定
    txt_u_32.config(state="normal")  # 一時的に編集可能に
    txt_u_32.delete(0, tk.END)
    txt_u_32.insert(0, round(UBtable['U'][2,1],4))
    txt_u_32.config(state="readonly") # 編集不可に設定
    txt_u_33.config(state="normal")  # 一時的に編集可能に
    txt_u_33.delete(0, tk.END)
    txt_u_33.insert(0, round(UBtable['U'][2,2],4))
    txt_u_33.config(state="readonly") # 編集不可に設定
    
    # B matrixを表示
    txt_b_11.config(state="normal")  # 一時的に編集可能に
    txt_b_11.delete(0, tk.END)
    txt_b_11.insert(0, round(UBtable['B'][0,0],4))
    txt_b_11.config(state="readonly") # 編集不可に設定
    txt_b_12.config(state="normal")  # 一時的に編集可能に
    txt_b_12.delete(0, tk.END)
    txt_b_12.insert(0, round(UBtable['B'][0,1],4))
    txt_b_12.config(state="readonly") # 編集不可に設定
    txt_b_13.config(state="normal")  # 一時的に編集可能に
    txt_b_13.delete(0, tk.END)
    txt_b_13.insert(0, round(UBtable['B'][0,2],4))
    txt_b_13.config(state="readonly") # 編集不可に設定
    txt_b_21.config(state="normal")  # 一時的に編集可能に
    txt_b_21.delete(0, tk.END)
    txt_b_21.insert(0, round(UBtable['B'][1,0],4))
    txt_b_21.config(state="readonly") # 編集不可に設定
    txt_b_22.config(state="normal")  # 一時的に編集可能に
    txt_b_22.delete(0, tk.END)
    txt_b_22.insert(0, round(UBtable['B'][1,1],4))
    txt_b_22.config(state="readonly") # 編集不可に設定
    txt_b_23.config(state="normal")  # 一時的に編集可能に
    txt_b_23.delete(0, tk.END)
    txt_b_23.insert(0, round(UBtable['B'][1,2],4))
    txt_b_23.config(state="readonly") # 編集不可に設定
    txt_b_31.config(state="normal")  # 一時的に編集可能に
    txt_b_31.delete(0, tk.END)
    txt_b_31.insert(0, round(UBtable['B'][2,0],4))
    txt_b_31.config(state="readonly") # 編集不可に設定
    txt_b_32.config(state="normal")  # 一時的に編集可能に
    txt_b_32.delete(0, tk.END)
    txt_b_32.insert(0, round(UBtable['B'][2,1],4))
    txt_b_32.config(state="readonly") # 編集不可に設定
    txt_b_33.config(state="normal")  # 一時的に編集可能に
    txt_b_33.delete(0, tk.END)
    txt_b_33.insert(0, round(UBtable['B'][2,2],4))
    txt_b_33.config(state="readonly") # 編集不可に設定
    
    # UB matrixを表示
    txt_ub_11.config(state="normal")  # 一時的に編集可能に
    txt_ub_11.delete(0, tk.END)
    txt_ub_11.insert(0, round(UBtable['UB'][0,0],4))
    txt_ub_11.config(state="readonly") # 編集不可に設定
    txt_ub_12.config(state="normal")  # 一時的に編集可能に
    txt_ub_12.delete(0, tk.END)
    txt_ub_12.insert(0, round(UBtable['UB'][0,1],4))
    txt_ub_12.config(state="readonly") # 編集不可に設定
    txt_ub_13.config(state="normal")  # 一時的に編集可能に
    txt_ub_13.delete(0, tk.END)
    txt_ub_13.insert(0, round(UBtable['UB'][0,2],4))
    txt_ub_13.config(state="readonly") # 編集不可に設定
    txt_ub_21.config(state="normal")  # 一時的に編集可能に
    txt_ub_21.delete(0, tk.END)
    txt_ub_21.insert(0, round(UBtable['UB'][1,0],4))
    txt_ub_21.config(state="readonly") # 編集不可に設定
    txt_ub_22.config(state="normal")  # 一時的に編集可能に
    txt_ub_22.delete(0, tk.END)
    txt_ub_22.insert(0, round(UBtable['UB'][1,1],4))
    txt_ub_22.config(state="readonly") # 編集不可に設定
    txt_ub_23.config(state="normal")  # 一時的に編集可能に
    txt_ub_23.delete(0, tk.END)
    txt_ub_23.insert(0, round(UBtable['UB'][1,2],4))
    txt_ub_23.config(state="readonly") # 編集不可に設定
    txt_ub_31.config(state="normal")  # 一時的に編集可能に
    txt_ub_31.delete(0, tk.END)
    txt_ub_31.insert(0, round(UBtable['UB'][2,0],4))
    txt_ub_31.config(state="readonly") # 編集不可に設定
    txt_ub_32.config(state="normal")  # 一時的に編集可能に
    txt_ub_32.delete(0, tk.END)
    txt_ub_32.insert(0, round(UBtable['UB'][2,1],4))
    txt_ub_32.config(state="readonly") # 編集不可に設定
    txt_ub_33.config(state="normal")  # 一時的に編集可能に
    txt_ub_33.delete(0, tk.END)
    txt_ub_33.insert(0, round(UBtable['UB'][2,2],4))
    txt_ub_33.config(state="readonly") # 編集不可に設定
    
    # Bragg peak positionの取得
    bpe = float(Energy.get())
    bpc2 = float(bp_c2.get())
    bpmu = float(bp_mu.get())
    bpnu = float(bp_nu.get())
    bp = np.array([float(bp_h.get()), float(bp_k.get()), float(bp_l.get())])
    
    # calculation pointの取得
    hw_ini = float(cqsef.get())
    hw_fin = float(cqset.get())
    hw_inc = float(cqsei.get())
    #h_cal = float(cqse1.get())
    #k_cal = float(cqse2.get())
    #l_cal = float(cqse3.get())
    
    # 分数も少数も受け付ける。
    # エントリーボックスから値を取得
    h_cal_entry = cqse1.get()
    try:
        # 少数の場合
        if '/' not in h_cal_entry:
            h_cal = float(h_cal_entry)
        else:
            # 分数の場合
            h_cal = float(Fraction(h_cal_entry))
    except ValueError:
        pass
    
    k_cal_entry = cqse2.get()
    try:
        # 少数の場合
        if '/' not in k_cal_entry:
            k_cal = float(k_cal_entry)
        else:
            # 分数の場合
            k_cal = float(Fraction(k_cal_entry))
    except ValueError:
        pass
    
    l_cal_entry = cqse3.get()
    try:
        # 少数の場合
        if '/' not in l_cal_entry:
            l_cal = float(l_cal_entry)
        else:
            # 分数の場合
            l_cal = float(Fraction(l_cal_entry))
    except ValueError:
        pass
               
    # Ei or Ef fixの判定
    fixe=float(eief.get())
    
    global angletable2
    angletable2 = angle_calc2(astar, bstar, cstar, U,B,UB, bpe, bpc2, bpmu, bpnu, bp, fixe, hw_ini, hw_fin, hw_inc, h_cal, k_cal, l_cal)
    
    # 先にフロアの条件を読み込み
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
    
    monochromator_radius = float(config['settings']['monochromator_radius'])
    monochromator_to_sample = float(config['settings']['monochromator_to_sample'])
    sample_goniometer_radius = float(config['settings']['sample_goniometer_radius'])
    sample_to_analyzer = float(config['settings']['sample_to_analyzer'])
    analyzer_radius = float(config['settings']['analyzer_radius'])
    analyzer_to_detector = float(config['settings']['analyzer_to_detector'])
    detector_radius = float(config['settings']['detector_radius'])
    floor_length = float(config['settings']['floor_length'])
    floor_width = float(config['settings']['floor_width'])
    floor_position_x = float(config['settings']['floor_position_x'])
    floor_position_y = float(config['settings']['floor_position_y'])
    
    A_sets = []  # A_setsリストを初期化
    QE_sets = []
    # resultsリストの各結果をTreeviewに追加
    for results in angletable2:
        values = tuple(results.values())
        item_id = tree.insert("", "end", values=values)
        # A1, A2, A3 を取得して A_sets に追加
        A1 = round(results['A1'], 4)  # 'A1'
        A2 = -round(results['A2'], 4)  # 'A2'
        A3 = round(results['A3'], 4)  # 'A3'
        A_sets.append([A1, A2, A3])  # A_sets に追加
        # hw, h,k,l
        hw = round(results['hw'], 4)  # 'A1'
        h = round(results['h'], 4)  # 'A2'
        k = round(results['k'], 4)  # 'A3'
        l = round(results['l'], 4)  # 'A3'
        QE_sets.append([hw, h, k,l])
        
        # sample gonioがfloorからはみ出る場合
        positionY_sample = monochromator_to_sample * np.sin(np.radians(results['A1'])) - sample_goniometer_radius # < floor_position_y
        # analyzer dramがfloorからはみ出る場合
        positionY_analyzer = monochromator_to_sample * np.sin(np.radians(results['A1'])) - sample_to_analyzer * np.sin(np.radians(results['A2'] - results['A1'])) - analyzer_radius # < floor_position_y
        # detector dramがfloorからはみ出る場合
        positionY_detector = monochromator_to_sample * np.sin(np.radians(results['A1'])) - sample_to_analyzer * np.sin(np.radians(results['A2'] - results['A1'])) + analyzer_to_detector * np.sin(np.radians(results['A3'] - results['A2'] + results['A1'])) - detector_radius # < floor_position_y
        positionX_detector = monochromator_to_sample * np.cos(np.radians(results['A1'])) + sample_to_analyzer * np.cos(np.radians(results['A2'] - results['A1'])) + analyzer_to_detector * np.cos(np.radians(results['A3'] - results['A2'] + results['A1'])) - detector_radius # < floor_position_x+floor_length
        
        # フロアからはみ出た場合、行を赤色にする
        if (positionY_sample < floor_position_y or
          positionY_analyzer < floor_position_y or
          positionY_detector < floor_position_y or
          positionX_detector > floor_position_x+floor_length):
            tree.tag_configure("red", foreground="red")  # 'red' タグを設定
            tree.item(item_id, tags=("red",))  # 行に 'red' タグを適用
        
        # hardware limitを超えた場合、行をい色にする
        elif (round(results['C1'],4)<float(hwl2f.get()) or
            round(results['C1'],4)>float(hwl2t.get()) or
            round(results['A1'],4)<float(hwl3f.get()) or
            round(results['A1'],4)>float(hwl3t.get()) or
            round(results['C2'],4)<float(hwl4f.get()) or
            round(results['C2'],4)>float(hwl4t.get()) or
            round(results['A2'],4)<float(hwl5f.get()) or
            round(results['A2'],4)>float(hwl5t.get()) or
            round(results['C3'],4)<float(hwl6f.get()) or
            round(results['C3'],4)>float(hwl6t.get()) or
            round(results['A3'],4)<float(hwl7f.get()) or
            round(results['A3'],4)>float(hwl7t.get()) or
            round(results['mu'],4)<float(hwl8f.get()) or
            round(results['mu'],4)>float(hwl8t.get()) or
            round(results['nu'],4)<float(hwl9f.get()) or
            round(results['nu'],4)>float(hwl9t.get())):
            tree.tag_configure("blue", foreground="blue")  # 'red' タグを設定
            tree.item(item_id, tags=("blue",))  # 行に 'red' タグを適用
        #tree.insert("", "end", values=values)
        
    if fig_spec.get()==1:
        # プロット関数を呼び出し
        #plot_spectrometer_with_gif(A_sets,QE_sets)
        plot_spectrometer(A_sets,QE_sets)

    Hfocus = calc_hf.get()
    num_ana = float(acna.get())
    
    if fig_reso.get()==1:
        bpe = float(Energy.get())
        fixe=float(eief.get())
        
        # Entry ウィジェットの値を辞書にまとめる
        entry_values = {
            "div_1st_m": div_1st_m.get(),
            "div_2nd_h": div_2nd_h.get(),
            "div_2nd_v": div_2nd_v.get(),
            "div_3rd_h": div_3rd_h.get(),
            "div_3rd_v": div_3rd_v.get(),
            "div_4th_h": div_4th_h.get(),
            "div_4th_v": div_4th_v.get(),
            "mos_mono_h": mos_mono_h.get(),
            "mos_mono_v": mos_mono_v.get(),
            "mos_sam_h": mos_sam_h.get(),
            "mos_sam_v": mos_sam_v.get(),
            "mos_ana_h": mos_ana_h.get(),
            "mos_ana_v": mos_ana_v.get(),
        }
        calcresolution_scan(A_sets,QE_sets,bpe,fixe,Hfocus,num_ana,entry_values)

    plt.show()
    
    return angletable2

# ボタンの作成
button = ttk.Button(tab_002a, text="calc", command=constQscan_show_table)
button.grid(row=4, column=0,columnspan=5, sticky="NSEW")

tab_002b = ttk.Labelframe(tab_002,text= "constant E scan")
tab_002b.grid(row=0,column=1,sticky="NSEW")
tab_002b.columnconfigure(0, weight=1)
tab_002b.columnconfigure(1, weight=1)
tab_002b.columnconfigure(2, weight=1)
tab_002b.columnconfigure(3, weight=1)
tab_002b.columnconfigure(4, weight=1)
tab_002b.rowconfigure(0, weight=1)
tab_002b.rowconfigure(1, weight=1)
tab_002b.rowconfigure(2, weight=1)
tab_002b.rowconfigure(3, weight=1)
tab_002b.rowconfigure(4, weight=1)

cesel= tk.Label(tab_002b,text='ℏω')
cesel.grid(row=0, column=4,sticky="NSEW")

cesl1= tk.Label(tab_002b,text='from',width=15)
cesl1.grid(row=1, column=0,sticky="NSEW")
cesl2= tk.Label(tab_002b,text='to')
cesl2.grid(row=2, column=0,sticky="NSEW")
cesl3= tk.Label(tab_002b,text='inc')
cesl3.grid(row=3, column=0,sticky="NSEW")

ceshl= tk.Label(tab_002b,text='h')
ceshl.grid(row=0, column=1,sticky="NSEW")
ceskl= tk.Label(tab_002b,text='k')
ceskl.grid(row=0, column=2,sticky="NSEW")
cesll= tk.Label(tab_002b,text='l')
cesll.grid(row=0, column=3,sticky="NSEW")

ces1 = ttk.Entry(tab_002b)
ces1.grid(row=1, column=1,sticky="NSEW")
ces1.insert(0,'-1')
ces2 = ttk.Entry(tab_002b)
ces2.grid(row=1, column=2,sticky="NSEW")
ces2.insert(0,'2')
ces3 = ttk.Entry(tab_002b)
ces3.grid(row=1, column=3,sticky="NSEW")
ces3.insert(0,'0')

ces4 = ttk.Entry(tab_002b)
ces4.grid(row=2, column=1,sticky="NSEW")
ces4.insert(0,'3')
ces5 = ttk.Entry(tab_002b)
ces5.grid(row=2, column=2,sticky="NSEW")
ces5.insert(0,'0')
ces6 = ttk.Entry(tab_002b)
ces6.grid(row=2, column=3,sticky="NSEW")
ces6.insert(0,'0')

ces7 = ttk.Entry(tab_002b)
ces7.grid(row=3, column=1,sticky="NSEW")
ces7.insert(0,'0.2')
ces8 = ttk.Entry(tab_002b)
ces8.grid(row=3, column=2,sticky="NSEW")
ces8.insert(0,'-0.1')
ces9 = ttk.Entry(tab_002b)
ces9.grid(row=3, column=3,sticky="NSEW")
ces9.insert(0,'0')

ces10 = ttk.Entry(tab_002b)
ces10.grid(row=1, column=4,sticky="NSEW")
ces10.insert(0,'1')

def conostEscan_show_table():
    # 新しいウィンドウを作成
    result_window = tk.Toplevel()
    result_window.title("calculation results (unit : deg)")
    
    # Treeviewの設定
    tree = ttk.Treeview(result_window, columns=("hw","h","k","l","C1", "A1", "C2", "A2", "C3", "A3", "mu", "nu"), show="headings")
    tree.pack(fill="both", expand=True)
    
    # 各列に見出しを設定
    for col in tree["columns"]:
        tree.heading(col, text=col)
        tree.column(col, width=80, anchor="center")
    
    # RLtableを取得し、辞書から必要な変数を取り出す
    RLtable = on_Rlcalc()
    astar = RLtable['astar']
    bstar = RLtable['bstar']
    cstar = RLtable['cstar']
    
    # UBtableを取得し、辞書から必要な変数を取り出す
    UBtable = on_UBcalc()
    
    U=UBtable['U']
    B=UBtable['B']
    UB=UBtable['UB']
    
    # displayに表示
    # U matrixを表示
    txt_u_11.config(state="normal")  # 一時的に編集可能に
    txt_u_11.delete(0, tk.END)
    txt_u_11.insert(0, round(UBtable['U'][0,0],4))
    txt_u_11.config(state="readonly") # 編集不可に設定
    txt_u_12.config(state="normal")  # 一時的に編集可能に
    txt_u_12.delete(0, tk.END)
    txt_u_12.insert(0, round(UBtable['U'][0,1],4))
    txt_u_12.config(state="readonly") # 編集不可に設定
    txt_u_13.config(state="normal")  # 一時的に編集可能に
    txt_u_13.delete(0, tk.END)
    txt_u_13.insert(0, round(UBtable['U'][0,2],4))
    txt_u_13.config(state="readonly") # 編集不可に設定
    txt_u_21.config(state="normal")  # 一時的に編集可能に
    txt_u_21.delete(0, tk.END)
    txt_u_21.insert(0, round(UBtable['U'][1,0],4))
    txt_u_21.config(state="readonly") # 編集不可に設定
    txt_u_22.config(state="normal")  # 一時的に編集可能に
    txt_u_22.delete(0, tk.END)
    txt_u_22.insert(0, round(UBtable['U'][1,1],4))
    txt_u_22.config(state="readonly") # 編集不可に設定
    txt_u_23.config(state="normal")  # 一時的に編集可能に
    txt_u_23.delete(0, tk.END)
    txt_u_23.insert(0, round(UBtable['U'][1,2],4))
    txt_u_23.config(state="readonly") # 編集不可に設定
    txt_u_31.config(state="normal")  # 一時的に編集可能に
    txt_u_31.delete(0, tk.END)
    txt_u_31.insert(0, round(UBtable['U'][2,0],4))
    txt_u_31.config(state="readonly") # 編集不可に設定
    txt_u_32.config(state="normal")  # 一時的に編集可能に
    txt_u_32.delete(0, tk.END)
    txt_u_32.insert(0, round(UBtable['U'][2,1],4))
    txt_u_32.config(state="readonly") # 編集不可に設定
    txt_u_33.config(state="normal")  # 一時的に編集可能に
    txt_u_33.delete(0, tk.END)
    txt_u_33.insert(0, round(UBtable['U'][2,2],4))
    txt_u_33.config(state="readonly") # 編集不可に設定
    
    # B matrixを表示
    txt_b_11.config(state="normal")  # 一時的に編集可能に
    txt_b_11.delete(0, tk.END)
    txt_b_11.insert(0, round(UBtable['B'][0,0],4))
    txt_b_11.config(state="readonly") # 編集不可に設定
    txt_b_12.config(state="normal")  # 一時的に編集可能に
    txt_b_12.delete(0, tk.END)
    txt_b_12.insert(0, round(UBtable['B'][0,1],4))
    txt_b_12.config(state="readonly") # 編集不可に設定
    txt_b_13.config(state="normal")  # 一時的に編集可能に
    txt_b_13.delete(0, tk.END)
    txt_b_13.insert(0, round(UBtable['B'][0,2],4))
    txt_b_13.config(state="readonly") # 編集不可に設定
    txt_b_21.config(state="normal")  # 一時的に編集可能に
    txt_b_21.delete(0, tk.END)
    txt_b_21.insert(0, round(UBtable['B'][1,0],4))
    txt_b_21.config(state="readonly") # 編集不可に設定
    txt_b_22.config(state="normal")  # 一時的に編集可能に
    txt_b_22.delete(0, tk.END)
    txt_b_22.insert(0, round(UBtable['B'][1,1],4))
    txt_b_22.config(state="readonly") # 編集不可に設定
    txt_b_23.config(state="normal")  # 一時的に編集可能に
    txt_b_23.delete(0, tk.END)
    txt_b_23.insert(0, round(UBtable['B'][1,2],4))
    txt_b_23.config(state="readonly") # 編集不可に設定
    txt_b_31.config(state="normal")  # 一時的に編集可能に
    txt_b_31.delete(0, tk.END)
    txt_b_31.insert(0, round(UBtable['B'][2,0],4))
    txt_b_31.config(state="readonly") # 編集不可に設定
    txt_b_32.config(state="normal")  # 一時的に編集可能に
    txt_b_32.delete(0, tk.END)
    txt_b_32.insert(0, round(UBtable['B'][2,1],4))
    txt_b_32.config(state="readonly") # 編集不可に設定
    txt_b_33.config(state="normal")  # 一時的に編集可能に
    txt_b_33.delete(0, tk.END)
    txt_b_33.insert(0, round(UBtable['B'][2,2],4))
    txt_b_33.config(state="readonly") # 編集不可に設定
    
    # UB matrixを表示
    txt_ub_11.config(state="normal")  # 一時的に編集可能に
    txt_ub_11.delete(0, tk.END)
    txt_ub_11.insert(0, round(UBtable['UB'][0,0],4))
    txt_ub_11.config(state="readonly") # 編集不可に設定
    txt_ub_12.config(state="normal")  # 一時的に編集可能に
    txt_ub_12.delete(0, tk.END)
    txt_ub_12.insert(0, round(UBtable['UB'][0,1],4))
    txt_ub_12.config(state="readonly") # 編集不可に設定
    txt_ub_13.config(state="normal")  # 一時的に編集可能に
    txt_ub_13.delete(0, tk.END)
    txt_ub_13.insert(0, round(UBtable['UB'][0,2],4))
    txt_ub_13.config(state="readonly") # 編集不可に設定
    txt_ub_21.config(state="normal")  # 一時的に編集可能に
    txt_ub_21.delete(0, tk.END)
    txt_ub_21.insert(0, round(UBtable['UB'][1,0],4))
    txt_ub_21.config(state="readonly") # 編集不可に設定
    txt_ub_22.config(state="normal")  # 一時的に編集可能に
    txt_ub_22.delete(0, tk.END)
    txt_ub_22.insert(0, round(UBtable['UB'][1,1],4))
    txt_ub_22.config(state="readonly") # 編集不可に設定
    txt_ub_23.config(state="normal")  # 一時的に編集可能に
    txt_ub_23.delete(0, tk.END)
    txt_ub_23.insert(0, round(UBtable['UB'][1,2],4))
    txt_ub_23.config(state="readonly") # 編集不可に設定
    txt_ub_31.config(state="normal")  # 一時的に編集可能に
    txt_ub_31.delete(0, tk.END)
    txt_ub_31.insert(0, round(UBtable['UB'][2,0],4))
    txt_ub_31.config(state="readonly") # 編集不可に設定
    txt_ub_32.config(state="normal")  # 一時的に編集可能に
    txt_ub_32.delete(0, tk.END)
    txt_ub_32.insert(0, round(UBtable['UB'][2,1],4))
    txt_ub_32.config(state="readonly") # 編集不可に設定
    txt_ub_33.config(state="normal")  # 一時的に編集可能に
    txt_ub_33.delete(0, tk.END)
    txt_ub_33.insert(0, round(UBtable['UB'][2,2],4))
    txt_ub_33.config(state="readonly") # 編集不可に設定
    
    # Bragg peak positionの取得
    bpe = float(Energy.get())
    bpc2 = float(bp_c2.get())
    bpmu = float(bp_mu.get())
    bpnu = float(bp_nu.get())
    bp = np.array([float(bp_h.get()), float(bp_k.get()), float(bp_l.get())])
    
    # calculation pointの取得
    hw_cal = float(ces10.get())
    #h_ini = float(ces1.get())
    #k_ini = float(ces2.get())
    #l_ini = float(ces3.get())
    #h_fin = float(ces4.get())
    #k_fin = float(ces5.get())
    #l_fin = float(ces6.get())
    #h_inc = float(ces7.get())
    #k_inc = float(ces8.get())
    #l_inc = float(ces9.get())
    
    h_ini_entry = ces1.get()
    try:
        # 少数の場合
        if '/' not in h_ini_entry:
            h_ini = float(h_ini_entry)
        else:
            # 分数の場合
            h_ini = float(Fraction(h_ini_entry))
    except ValueError:
        pass
    
    k_ini_entry = ces2.get()
    try:
        # 少数の場合
        if '/' not in k_ini_entry:
            k_ini = float(k_ini_entry)
        else:
            # 分数の場合
            k_ini = float(Fraction(k_ini_entry))
    except ValueError:
        pass
    
    l_ini_entry = ces3.get()
    try:
        # 少数の場合
        if '/' not in l_ini_entry:
            l_ini = float(l_ini_entry)
        else:
            # 分数の場合
            l_ini = float(Fraction(l_ini_entry))
    except ValueError:
        pass
    
    h_fin_entry = ces4.get()
    try:
        # 少数の場合
        if '/' not in h_fin_entry:
            h_fin = float(h_fin_entry)
        else:
            # 分数の場合
            h_fin = float(Fraction(h_fin_entry))
    except ValueError:
        pass
    
    k_fin_entry = ces5.get()
    try:
        # 少数の場合
        if '/' not in k_fin_entry:
            k_fin = float(k_fin_entry)
        else:
            # 分数の場合
            k_fin = float(Fraction(k_fin_entry))
    except ValueError:
        pass
    
    l_fin_entry = ces6.get()
    try:
        # 少数の場合
        if '/' not in l_fin_entry:
            l_fin = float(l_fin_entry)
        else:
            # 分数の場合
            l_fin = float(Fraction(l_fin_entry))
    except ValueError:
        pass
    
    h_inc_entry = ces7.get()
    try:
        # 少数の場合
        if '/' not in h_inc_entry:
            h_inc = float(h_inc_entry)
        else:
            # 分数の場合
            h_inc = float(Fraction(h_inc_entry))
    except ValueError:
        pass
    
    k_inc_entry = ces8.get()
    try:
        # 少数の場合
        if '/' not in k_inc_entry:
            k_inc = float(k_inc_entry)
        else:
            # 分数の場合
            k_inc = float(Fraction(k_inc_entry))
    except ValueError:
        pass
    
    l_inc_entry = ces9.get()
    try:
        # 少数の場合
        if '/' not in l_inc_entry:
            l_inc = float(l_inc_entry)
        else:
            # 分数の場合
            l_inc = float(Fraction(l_inc_entry))
    except ValueError:
        pass
    
    # 先にフロアの条件を読み込み
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
    
    monochromator_radius = float(config['settings']['monochromator_radius'])
    monochromator_to_sample = float(config['settings']['monochromator_to_sample'])
    sample_goniometer_radius = float(config['settings']['sample_goniometer_radius'])
    sample_to_analyzer = float(config['settings']['sample_to_analyzer'])
    analyzer_radius = float(config['settings']['analyzer_radius'])
    analyzer_to_detector = float(config['settings']['analyzer_to_detector'])
    detector_radius = float(config['settings']['detector_radius'])
    floor_length = float(config['settings']['floor_length'])
    floor_width = float(config['settings']['floor_width'])
    floor_position_x = float(config['settings']['floor_position_x'])
    floor_position_y = float(config['settings']['floor_position_y'])

    # Ei or Ef fixの判定
    fixe=float(eief.get())
    
    global angletable3
    angletable3 = angle_calc3(astar,bstar,cstar,U,B,UB,bpe,bpc2,bpmu,bpnu,bp,fixe,hw_cal,h_ini,k_ini,l_ini,h_fin,k_fin,l_fin,h_inc,k_inc,l_inc)
    
    A_sets = []  # A_setsリストを初期化
    QE_sets = []
    # resultsリストの各結果をTreeviewに追加
    for results in angletable3:
        values = tuple(results.values())
        item_id = tree.insert("", "end", values=values)
    
        # A1, A2, A3 を取得して A_sets に追加
        A1 = round(results['A1'], 4)  # 'A1'
        A2 = -round(results['A2'], 4)  # 'A2'
        A3 = round(results['A3'], 4)  # 'A3'
        A_sets.append([A1, A2, A3])  # A_sets に追加
        # hw, h,k,l
        hw = round(results['hw'], 4)  # 'A1'
        h = round(results['h'], 4)  # 'A2'
        k = round(results['k'], 4)  # 'A3'
        l = round(results['l'], 4)  # 'A3'
        QE_sets.append([hw, h, k,l])
        
        # sample gonioがfloorからはみ出る場合
        positionY_sample = monochromator_to_sample * np.sin(np.radians(results['A1'])) - sample_goniometer_radius # < floor_position_y
        # analyzer dramがfloorからはみ出る場合
        positionY_analyzer = monochromator_to_sample * np.sin(np.radians(results['A1'])) - sample_to_analyzer * np.sin(np.radians(results['A2'] - results['A1'])) - analyzer_radius # < floor_position_y
        # detector dramがfloorからはみ出る場合
        positionY_detector = monochromator_to_sample * np.sin(np.radians(results['A1'])) - sample_to_analyzer * np.sin(np.radians(results['A2'] - results['A1'])) + analyzer_to_detector * np.sin(np.radians(results['A3'] - results['A2'] + results['A1'])) - detector_radius # < floor_position_y
        positionX_detector = monochromator_to_sample * np.cos(np.radians(results['A1'])) + sample_to_analyzer * np.cos(np.radians(results['A2'] - results['A1'])) + analyzer_to_detector * np.cos(np.radians(results['A3'] - results['A2'] + results['A1'])) - detector_radius # < floor_position_x+floor_length
        
        # フロアからはみ出た場合、行を赤色にする
        if (positionY_sample < floor_position_y or
          positionY_analyzer < floor_position_y or
          positionY_detector < floor_position_y or
          positionX_detector > floor_position_x+floor_length):
            tree.tag_configure("red", foreground="red")  # 'red' タグを設定
            tree.item(item_id, tags=("red",))  # 行に 'red' タグを適用

        # hardware limitを超えた場合、行をい色にする
        elif (round(results['C1'],4)<float(hwl2f.get()) or
            round(results['C1'],4)>float(hwl2t.get()) or
            round(results['A1'],4)<float(hwl3f.get()) or
            round(results['A1'],4)>float(hwl3t.get()) or
            round(results['C2'],4)<float(hwl4f.get()) or
            round(results['C2'],4)>float(hwl4t.get()) or
            round(results['A2'],4)<float(hwl5f.get()) or
            round(results['A2'],4)>float(hwl5t.get()) or
            round(results['C3'],4)<float(hwl6f.get()) or
            round(results['C3'],4)>float(hwl6t.get()) or
            round(results['A3'],4)<float(hwl7f.get()) or
            round(results['A3'],4)>float(hwl7t.get()) or
            round(results['mu'],4)<float(hwl8f.get()) or
            round(results['mu'],4)>float(hwl8t.get()) or
            round(results['nu'],4)<float(hwl9f.get()) or
            round(results['nu'],4)>float(hwl9t.get())):
            tree.tag_configure("blue", foreground="blue")  # 'red' タグを設定
            tree.item(item_id, tags=("blue",))  # 行に 'red' タグを適用
        
    if fig_spec.get()==1:
        # プロット関数を呼び出し
        #plot_spectrometer_with_gif(A_sets,QE_sets)
        plot_spectrometer(A_sets,QE_sets)
    
    Hfocus = calc_hf.get()
    num_ana = float(acna.get())
    
    if fig_reso.get()==1:
        bpe = float(Energy.get())
        fixe=float(eief.get())
        
        # Entry ウィジェットの値を辞書にまとめる
        entry_values = {
            "div_1st_m": div_1st_m.get(),
            "div_2nd_h": div_2nd_h.get(),
            "div_2nd_v": div_2nd_v.get(),
            "div_3rd_h": div_3rd_h.get(),
            "div_3rd_v": div_3rd_v.get(),
            "div_4th_h": div_4th_h.get(),
            "div_4th_v": div_4th_v.get(),
            "mos_mono_h": mos_mono_h.get(),
            "mos_mono_v": mos_mono_v.get(),
            "mos_sam_h": mos_sam_h.get(),
            "mos_sam_v": mos_sam_v.get(),
            "mos_ana_h": mos_ana_h.get(),
            "mos_ana_v": mos_ana_v.get(),
        }
        calcresolution_scan(A_sets,QE_sets,bpe,fixe,Hfocus,num_ana,entry_values)
        
    plt.show()
    
    return angletable3

# ボタンの作成
button = ttk.Button(tab_002b, text="calc", command=conostEscan_show_table)
button.grid(row=4, column=0,columnspan=5, sticky="NSEW")

# グリッドの重みを設定
tab_003.columnconfigure(0, weight=1)
tab_003.rowconfigure(0, weight=2)
tab_003.rowconfigure(1, weight=3)

tab_003a = ttk.Labelframe(tab_003,text= "fitting lattice paramter")
tab_003a.grid(row=0,column=0,sticky="NSEW")
tab_003a.columnconfigure(0, weight=1)
tab_003a.columnconfigure(1, weight=1)
tab_003a.columnconfigure(2, weight=1)
tab_003a.columnconfigure(3, weight=1)
tab_003a.columnconfigure(4, weight=1)
tab_003a.columnconfigure(5, weight=1)
tab_003a.rowconfigure(0, weight=1)
tab_003a.rowconfigure(1, weight=1)

tab_003b = ttk.Labelframe(tab_003,text= "fitting results")
tab_003b.grid(row=1,column=0,sticky="NSEW")
tab_003b.columnconfigure(0, weight=1)
tab_003b.columnconfigure(1, weight=1)
tab_003b.columnconfigure(2, weight=1)
tab_003b.columnconfigure(3, weight=1)
tab_003b.columnconfigure(4, weight=1)
tab_003b.columnconfigure(5, weight=1)
tab_003b.columnconfigure(6, weight=1)
tab_003b.rowconfigure(0, weight=1)
tab_003b.rowconfigure(1, weight=1)
tab_003b.rowconfigure(2, weight=1)

# 結晶系コンボボックスのリスト
cryst_list =['cubic','tetragonal', 'orthorhombic', 'hexagonal' ,'monoclinic', 'triclinic']

# Xコンボボックスのラベル
lbl_cl = tk.Label(tab_003a,text='crystal')
lbl_cl.grid(row=0, column=0,sticky="NSEW")

# Xコンボボックスを設置
cb_cl = ttk.Combobox(tab_003a, values = cryst_list,width=18)
cb_cl.grid(row=1, column=0,sticky="NSEW")

#コンボボックスのリストの先頭を表示
cb_cl.set(cryst_list[0])

cry_fit= tk.Label(tab_003a,text='h')
cry_fit.grid(row=0, column=1,sticky="NSEW")
cry_fit_h = ttk.Entry(tab_003a)
cry_fit_h.grid(row=1, column=1,sticky="NSEW")
cry_fit_h.insert(0,'1.5')

cry_fit= tk.Label(tab_003a,text='k')
cry_fit.grid(row=0, column=2,sticky="NSEW")
cry_fit_k = ttk.Entry(tab_003a)
cry_fit_k.grid(row=1, column=2,sticky="NSEW")
cry_fit_k.insert(0,'0')

cry_fit= tk.Label(tab_003a,text='l')
cry_fit.grid(row=0, column=3,sticky="NSEW")
cry_fit_l = ttk.Entry(tab_003a)
cry_fit_l.grid(row=1, column=3,sticky="NSEW")
cry_fit_l.insert(0,'0')

cry_fit= tk.Label(tab_003a,text='A2')
cry_fit.grid(row=0, column=4,sticky="NSEW")
cry_fit_a2 = ttk.Entry(tab_003a)
cry_fit_a2.grid(row=1, column=4,sticky="NSEW")
cry_fit_a2.insert(0,'90.6')

def fitting_process():
    # selectされたindexを読み込む
    cry_select = cb_cl.current()
    
    # 変数初期化
    #h = float(cry_fit_h.get())     # h
    #k = float(cry_fit_k.get())     # k
    #l = float(cry_fit_l.get())     # l
    
    h_fit_entry = cry_fit_h.get()
    try:
        # 少数の場合
        if '/' not in h_fit_entry:
            h = float(h_fit_entry)
        else:
            # 分数の場合
            h = float(Fraction(h_fit_entry))
    except ValueError:
        pass
    
    k_fit_entry = cry_fit_k.get()
    try:
        # 少数の場合
        if '/' not in k_fit_entry:
            k = float(k_fit_entry)
        else:
            # 分数の場合
            k = float(Fraction(k_fit_entry))
    except ValueError:
        pass
    
    l_fit_entry = cry_fit_l.get()
    try:
        # 少数の場合
        if '/' not in l_fit_entry:
            l = float(l_fit_entry)
        else:
            # 分数の場合
            l = float(Fraction(l_fit_entry))
    except ValueError:
        pass
    
    kf = (float(Energy.get()) / 2.072) ** (1 / 2)
    ki = (float(Energy.get()) / 2.072) ** (1 / 2)
    
    # 結晶系の選択に基づいてフィッティング処理を分ける
    cryst_type = ''
    if cry_select == 0:  # cubic
        cryst_type = 'cubic'
    elif cry_select == 1:  # tetragonal
        cryst_type = 'tetragonal'
    elif cry_select == 2:  # orthorhombic
        cryst_type = 'orthorhombic'
    elif cry_select == 3:  # hexagonal
        cryst_type = 'hexagonal'
    elif cry_select == 4:  # monoclinic
        cryst_type = 'monoclinic'
    elif cry_select == 5:  # triclinic
        cryst_type = 'triclinic'
    
    # フィッティング処理を実行
    try:
        # a2_measuredは仮の値として設定していますが、実際のデータをここに渡してください
        a2_measured = float(cry_fit_a2.get())  # 実際の測定値を使う

        # GUI から初期値を取得
        initial_params = get_parameters()  # これが GUI の初期値を取得する部分
        
        # 修正済みの fit_lattice_constants を呼び出す
        result = fit_lattice_constants(
            ki, kf, (h, k, l), a2_measured, cryst_type,
            initial_params=initial_params  # 取得した初期値を渡す
        )
        
        # 結果をGUIに反映
        fit_resa.configure(state="normal")
        fit_resb.configure(state="normal")
        fit_resc.configure(state="normal")
        fit_resal.configure(state="normal")
        fit_resbe.configure(state="normal")
        fit_resga.configure(state="normal")

        fit_resa.delete(0, tk.END)
        fit_resb.delete(0, tk.END)
        fit_resc.delete(0, tk.END)
        fit_resal.delete(0, tk.END)
        fit_resbe.delete(0, tk.END)
        fit_resga.delete(0, tk.END)

        fit_resa.insert(0, round(result.get("a", "N/A"),4))
        fit_resb.insert(0, round(result.get("b", "N/A"),4))
        fit_resc.insert(0, round(result.get("c", "N/A"),4))
        fit_resal.insert(0, round(result.get("alpha", "N/A"),4))
        fit_resbe.insert(0, round(result.get("beta", "N/A"),4))
        fit_resga.insert(0, round(result.get("gamma", "N/A"),4))

        fit_resa.configure(state="readonly")
        fit_resb.configure(state="readonly")
        fit_resc.configure(state="readonly")
        fit_resal.configure(state="readonly")
        fit_resbe.configure(state="readonly")
        fit_resga.configure(state="readonly")

        # エラーメッセージラベルをクリア
        fit_resme.configure(text="The fitting was successful.", fg="black")

    except ValueError as e:
        # エラー発生時、エラーメッセージを表示
        fit_resme.configure(text=f"{str(e)}", fg="red")

# フィッティング開始ボタン
fit_button = ttk.Button(tab_003a, text="Fit", command=fitting_process,width=16)
fit_button.grid(row=1, column=5, sticky="NSEW")

# fitting結果の表示
fit_res1= tk.Label(tab_003b,text='a')
fit_res1.grid(row=0, column=0,sticky="NSEW")
fit_res2= tk.Label(tab_003b,text='b')
fit_res2.grid(row=0, column=1,sticky="NSEW")
fit_res3= tk.Label(tab_003b,text='c')
fit_res3.grid(row=0, column=2,sticky="NSEW")
fit_res4= tk.Label(tab_003b,text='α')
fit_res4.grid(row=0, column=3,sticky="NSEW")
fit_res5= tk.Label(tab_003b,text='β')
fit_res5.grid(row=0, column=4,sticky="NSEW")
fit_res6= tk.Label(tab_003b,text='γ')
fit_res6.grid(row=0, column=5,sticky="NSEW")

fit_resa = ttk.Entry(tab_003b,state="readonly")
fit_resa.grid(row=1, column=0,sticky="NSEW")
fit_resb = ttk.Entry(tab_003b,state="readonly")
fit_resb.grid(row=1, column=1,sticky="NSEW")
fit_resc = ttk.Entry(tab_003b,state="readonly")
fit_resc.grid(row=1, column=2,sticky="NSEW")
fit_resal = ttk.Entry(tab_003b,state="readonly")
fit_resal.grid(row=1, column=3,sticky="NSEW")
fit_resbe = ttk.Entry(tab_003b,state="readonly")
fit_resbe.grid(row=1, column=4,sticky="NSEW")
fit_resga = ttk.Entry(tab_003b,state="readonly")
fit_resga.grid(row=1, column=5,sticky="NSEW")

def reflection():
    fit_a =float(fit_resa.get())
    fit_b =float(fit_resb.get())
    fit_c =float(fit_resc.get())
    fit_al =float(fit_resal.get())
    fit_be =float(fit_resbe.get())
    fit_ga =float(fit_resga.get())
    
    la.delete(0, tk.END)
    lb.delete(0, tk.END)
    lc.delete(0, tk.END)
    lc_alpha.delete(0, tk.END)
    lc_beta.delete(0, tk.END)
    lc_gamma.delete(0, tk.END)
    
    la.insert(0, fit_a)
    lb.insert(0, fit_b)
    lc.insert(0, fit_c)
    lc_alpha.insert(0, fit_al)
    lc_beta.insert(0, fit_be)
    lc_gamma.insert(0, fit_ga)

# fitting結果を反映させるボタン
ref_button = ttk.Button(tab_003b, text="set param", command=reflection,width=18)
ref_button.grid(row=1, column=6, sticky="NSEW")

fit_reswa = tk.Label(tab_003b,text='warning : ')
fit_reswa.grid(row=2, column=0,sticky="NSEW")

fit_resme = tk.Label(tab_003b,text='')
fit_resme.grid(row=2, column=1,columnspan=5,sticky="NSEW")

# toolbox
# グリッドの重みを設定
tab_004.columnconfigure(0, weight=4)
tab_004.columnconfigure(1, weight=4)
tab_004.rowconfigure(0, weight=2)
tab_004.rowconfigure(1, weight=4)

tab_004a = ttk.Labelframe(tab_004,text= "reciprocal vectors")
tab_004a.grid(row=0,column=0,rowspan=2, sticky="NSEW")
tab_004a.columnconfigure(0, weight=1)
tab_004a.columnconfigure(1, weight=1)
tab_004a.columnconfigure(2, weight=1)
tab_004a.columnconfigure(3, weight=1)
tab_004a.rowconfigure(0, weight=1)
tab_004a.rowconfigure(1, weight=1)
tab_004a.rowconfigure(2, weight=1)
tab_004a.rowconfigure(3, weight=1)
tab_004a.rowconfigure(4, weight=1)

rvx= tk.Label(tab_004a,text='x')
rvx.grid(row=0, column=1,sticky="NSEW")
rvy= tk.Label(tab_004a,text='y')
rvy.grid(row=0, column=2,sticky="NSEW")
rvz= tk.Label(tab_004a,text='z')
rvz.grid(row=0, column=3,sticky="NSEW")

rv1= tk.Label(tab_004a,text='a* = ',width=16)
rv1.grid(row=1, column=0,sticky="NSEW")
rv1x = ttk.Entry(tab_004a,state="readonly")
rv1x.grid(row=1, column=1,sticky="NSEW")
rv1y = ttk.Entry(tab_004a,state="readonly")
rv1y.grid(row=1, column=2,sticky="NSEW")
rv1z = ttk.Entry(tab_004a,state="readonly")
rv1z.grid(row=1, column=3,sticky="NSEW")

rv2= tk.Label(tab_004a,text='b* = ',width=16)
rv2.grid(row=2, column=0,sticky="NSEW")
rv2x = ttk.Entry(tab_004a,state="readonly")
rv2x.grid(row=2, column=1,sticky="NSEW")
rv2y = ttk.Entry(tab_004a,state="readonly")
rv2y.grid(row=2, column=2,sticky="NSEW")
rv2z = ttk.Entry(tab_004a,state="readonly")
rv2z.grid(row=2, column=3,sticky="NSEW")

rv3= tk.Label(tab_004a,text='c* = ',width=16)
rv3.grid(row=3, column=0,sticky="NSEW")
rv3x = ttk.Entry(tab_004a,state="readonly")
rv3x.grid(row=3, column=1,sticky="NSEW")
rv3y = ttk.Entry(tab_004a,state="readonly")
rv3y.grid(row=3, column=2,sticky="NSEW")
rv3z = ttk.Entry(tab_004a,state="readonly")
rv3z.grid(row=3, column=3,sticky="NSEW")

def calc_RVs():
    # RLtableを取得し、辞書から必要な変数を取り出す
    RLtable = on_Rlcalc()
    
    # 結果をGUIに反映
    rv1x.configure(state="normal")
    rv1y.configure(state="normal")
    rv1z.configure(state="normal")
    rv1x.delete(0, tk.END)
    rv1y.delete(0, tk.END)
    rv1z.delete(0, tk.END)
    rv1x.insert(0, round(RLtable['astar'][0],4))
    rv1y.insert(0, round(RLtable['astar'][1],4))
    rv1z.insert(0, round(RLtable['astar'][2],4))
    rv1x.configure(state="readonly")
    rv1y.configure(state="readonly")
    rv1z.configure(state="readonly")
    
    rv2x.configure(state="normal")
    rv2y.configure(state="normal")
    rv2z.configure(state="normal")
    rv2x.delete(0, tk.END)
    rv2y.delete(0, tk.END)
    rv2z.delete(0, tk.END)
    rv2x.insert(0, round(RLtable['bstar'][0],4))
    rv2y.insert(0, round(RLtable['bstar'][1],4))
    rv2z.insert(0, round(RLtable['bstar'][2],4))
    rv2x.configure(state="readonly")
    rv2y.configure(state="readonly")
    rv2z.configure(state="readonly")
    
    rv3x.configure(state="normal")
    rv3y.configure(state="normal")
    rv3z.configure(state="normal")
    rv3x.delete(0, tk.END)
    rv3y.delete(0, tk.END)
    rv3z.delete(0, tk.END)
    rv3x.insert(0, round(RLtable['cstar'][0],4))
    rv3y.insert(0, round(RLtable['cstar'][1],4))
    rv3z.insert(0, round(RLtable['cstar'][2],4))
    rv3x.configure(state="readonly")
    rv3y.configure(state="readonly")
    rv3z.configure(state="readonly")

# fitting結果を反映させるボタン
ref_button = ttk.Button(tab_004a, text="calc", command=calc_RVs)
ref_button.grid(row=4, column=0, columnspan=4 , sticky="NSEW")

tab_004b = ttk.Labelframe(tab_004,text= "calculation1")
tab_004b.grid(row=0,column=1,sticky="NSEW")
tab_004b.columnconfigure(0, weight=1)
tab_004b.columnconfigure(1, weight=1)
tab_004b.columnconfigure(2, weight=1)
tab_004b.rowconfigure(0, weight=1)
tab_004b.rowconfigure(1, weight=1)

def trans_ELK(event=None):
    """
    入力された値に基づいてエネルギー、波長、波数を計算し、エントリーボックスを更新します。
    """
    try:
        # フォーカスされているウィジェットを特定
        focused_widget = root.focus_get()

        # 各エントリーボックスの値を取得
        energy = entry_energy.get().strip()
        wavelength = entry_wavelength.get().strip()
        wavenumber = entry_wavenumber.get().strip()

        # フォーカスされているエントリーボックスに応じて処理
        if focused_widget == entry_energy and energy:
            energy = float(energy)
            wavelength = math.sqrt(81.81 / energy)  # Å
            wavenumber = 2 * math.pi / wavelength  # Å⁻¹
            # 他のボックスをクリアして計算結果を出力
            entry_wavelength.delete(0, tk.END)
            entry_wavenumber.delete(0, tk.END)
            entry_wavelength.insert(0, f"{wavelength:.4f}")
            entry_wavenumber.insert(0, f"{wavenumber:.4f}")

        elif focused_widget == entry_wavelength and wavelength:
            wavelength = float(wavelength)
            energy = 81.81 / (wavelength ** 2)  # meV
            wavenumber = 2 * math.pi / wavelength  # Å⁻¹
            # 他のボックスをクリアして計算結果を出力
            entry_energy.delete(0, tk.END)
            entry_wavenumber.delete(0, tk.END)
            entry_energy.insert(0, f"{energy:.4f}")
            entry_wavenumber.insert(0, f"{wavenumber:.4f}")

        elif focused_widget == entry_wavenumber and wavenumber:
            wavenumber = float(wavenumber)
            wavelength = 2 * math.pi / wavenumber  # Å
            energy = 81.81 / (wavelength ** 2)  # meV
            # 他のボックスをクリアして計算結果を出力
            entry_energy.delete(0, tk.END)
            entry_wavelength.delete(0, tk.END)
            entry_wavelength.insert(0, f"{wavelength:.4f}")
            entry_energy.insert(0, f"{energy:.4f}")

        else:
            pass
            return

    except ValueError:
        pass

# ラベルとエントリーボックス
tk.Label(tab_004b, text="E (meV)").grid(row=0, column=0,sticky="NSEW")
entry_energy = ttk.Entry(tab_004b)
entry_energy.grid(row=1, column=0,sticky="NSEW")

tk.Label(tab_004b, text="λ (Å)").grid(row=0, column=1,sticky="NSEW")
entry_wavelength = ttk.Entry(tab_004b)
entry_wavelength.grid(row=1, column=1,sticky="NSEW")

tk.Label(tab_004b, text="k (Å⁻¹)").grid(row=0, column=2,sticky="NSEW")
entry_wavenumber = ttk.Entry(tab_004b)
entry_wavenumber.grid(row=1, column=2,sticky="NSEW")

# エンターキーで計算を実行
entry_energy.bind("<Return>", trans_ELK)
entry_wavelength.bind("<Return>", trans_ELK)
entry_wavenumber.bind("<Return>", trans_ELK)

tab_004c = ttk.Labelframe(tab_004,text= "calculation2")
tab_004c.grid(row=1,column=1,sticky="NSEW")
tab_004c.columnconfigure(0, weight=1)
tab_004c.columnconfigure(1, weight=1)
tab_004c.columnconfigure(2, weight=1)
tab_004c.columnconfigure(3, weight=1)
tab_004c.rowconfigure(0, weight=1)
tab_004c.rowconfigure(1, weight=1)
tab_004c.rowconfigure(2, weight=1)
tab_004c.rowconfigure(3, weight=1)

# ラベルとエントリーボックス
tk.Label(tab_004c, text="ℏω (meV)").grid(row=0, column=0,sticky="NSEW")
entry_hw = ttk.Entry(tab_004c)
entry_hw.grid(row=1, column=0,sticky="NSEW")

tk.Label(tab_004c, text="h").grid(row=0, column=1,sticky="NSEW")
entry_h = ttk.Entry(tab_004c)
entry_h.grid(row=1, column=1,sticky="NSEW")

tk.Label(tab_004c, text="k").grid(row=0, column=2,sticky="NSEW")
entry_k = ttk.Entry(tab_004c)
entry_k.grid(row=1, column=2,sticky="NSEW")

tk.Label(tab_004c, text="l").grid(row=0, column=3,sticky="NSEW")
entry_l = ttk.Entry(tab_004c)
entry_l.grid(row=1, column=3,sticky="NSEW")

# ttaとdとQを計算する
def calc_tta_d_Q():
    RLtable = on_Rlcalc()
    astar = RLtable['astar']
    bstar = RLtable['bstar']
    cstar = RLtable['cstar']
    
    entry_tta.configure(state="normal")
    entry_d.configure(state="normal")
    entry_q.configure(state="normal")
    entry_tta.delete(0,tk.END)
    entry_d.delete(0,tk.END)
    entry_q.delete(0,tk.END)
    
    hw=float(entry_hw.get())
    h=float(entry_h.get())
    k=float(entry_k.get())
    l=float(entry_l.get())
    hkl=h*astar+k*bstar+l*cstar

    #計算されたrlu
    Nhkl=np.linalg.norm(hkl)
    dhkl=2*math.pi/Nhkl
    
    # ei / ef
    if float(eief.get())==0: # ei fix
        Ei = float(Energy.get())
        Ef = float(Energy.get()) - hw
    elif float(eief.get())==1: # ef fix
        Ei = float(Energy.get()) + hw
        Ef = float(Energy.get())

    #lamdaとkを計算
    Ki=(Ei/2.072)**(1/2)
    Kf=(Ef/2.072)**(1/2)
    
    tta=math.degrees(math.acos((Ki**2+Kf**2-Nhkl**2)/(2*Ki*Kf)))
    
    entry_tta.insert(0, f"{tta:.4f}")
    entry_d.insert(0, f"{dhkl:.4f}")
    entry_q.insert(0, f"{Nhkl:.4f}")
    
    entry_tta.configure(state="readonly")
    entry_d.configure(state="readonly")
    entry_q.configure(state="readonly")
    
# 結果を反映させるボタン
ref_button = ttk.Button(tab_004c, text="calc", command=calc_tta_d_Q)
ref_button.grid(row=3, column=3 , sticky="NSEW")

# ラベルとエントリーボックス
tk.Label(tab_004c, text="2θ (deg)").grid(row=2, column=0,sticky="NSEW")
entry_tta = ttk.Entry(tab_004c,state="readonly")
entry_tta.grid(row=3, column=0,sticky="NSEW")

tk.Label(tab_004c, text="d (Å)").grid(row=2, column=1,sticky="NSEW")
entry_d = ttk.Entry(tab_004c,state="readonly")
entry_d.grid(row=3, column=1,sticky="NSEW")

tk.Label(tab_004c, text="Q (Å⁻¹)").grid(row=2, column=2,sticky="NSEW")
entry_q = ttk.Entry(tab_004c,state="readonly")
entry_q.grid(row=3, column=2,sticky="NSEW")

#メニューバーの作成
menubar = tk.Menu(root)
root.configure(menu=menubar)

#fileメニュー(setting)
filemenu = tk.Menu(menubar,tearoff=0)
menubar.add_cascade(label="ini.file",menu=filemenu)
#fileメニューにini fileのload
filemenu.add_command(label="load ini.file",command=load_values_from_ini)
#fileメニューにexitを追加。ついでにexit funcも実装

#fileメニューにini fileのsave
filemenu.add_command(label="save ini.file",command=save_values_to_ini)

def save_cQ_table():

    # 保存ダイアログを表示してファイル名を取得
    file_path = filedialog.asksaveasfilename(
        defaultextension=".csv",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        title="Save as CSV"
    )

    if file_path:
        with open(file_path, mode='w', newline='') as file:
            writer = csv.writer(file)

            # ヘッダーを書き込む
            header = ['hw','h','k','l','C1', 'A1', 'C2', 'A2', 'C3', 'A3', 'mu', 'nu']  # ヘッダー名を必要に応じて調整
            writer.writerow(header)

            # angletable2 の各結果を CSV に書き込む
            for results in angletable2:
                # results は辞書なので、values() で値だけ取り出してタプルにする
                values = tuple(results.values())
                writer.writerow(values)

def save_cE_table():

    # 保存ダイアログを表示してファイル名を取得
    file_path = filedialog.asksaveasfilename(
        defaultextension=".csv",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        title="Save as CSV"
    )

    if file_path:
        with open(file_path, mode='w', newline='') as file:
            writer = csv.writer(file)

            # ヘッダーを書き込む
            header = ['hw','h','k','l','C1', 'A1', 'C2', 'A2', 'C3', 'A3', 'mu', 'nu']  # ヘッダー名を必要に応じて調整
            writer.writerow(header)

            # angletable2 の各結果を CSV に書き込む
            for results in angletable3:
                # results は辞書なので、values() で値だけ取り出してタプルにする
                values = tuple(results.values())
                writer.writerow(values)

#fileメニュー(save)
filemenu2 = tk.Menu(menubar,tearoff=0)
menubar.add_cascade(label="output",menu=filemenu2)
#fileメニュー(save)
filemenu2.add_command(label="const Q scan",command=save_cQ_table)
#fileメニュー(save)
filemenu2.add_command(label="const E scan",command=save_cE_table)

#fileメニュー(exit)
filemenu3 = tk.Menu(menubar,tearoff=0)
menubar.add_cascade(label="exit",menu=filemenu3)
filemenu3.add_command(label="exit",command=lambda:root.destroy())

# アプリ起動時にデフォルト値を読み込む
load_values_from_ini()

#window状態の維持
root.mainloop()

#############
# pyinstaller tips
# cd C:\DATA_HK\python\SPICE_python
# pyinstaller --noconsole --onefile --add-data "config.ini;." --add-data "logo2.ico;." --icon=logo2.ico TriAxionSim.py

"""
コマンド
pyinstaller TriAxionSim.py --noconsole
--onedir or -D
出力を1ディレクトリにまとめる
--onefile or -F
出力を1ファイルにまとめる
--noconsole or -w
コンソールを表示しない
--clean
ビルド前に前回のキャッシュと出力ディレクトリを削除
"""
#############
# %%
