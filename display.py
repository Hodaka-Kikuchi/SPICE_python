# Automation System for Utility in Research Analysis
# ASYURA
# UB matrixによるoffsetの自動入力
# 右上にバージョン情報を表示
__version__ = '1.0.0'
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

import statistics as stat
from scipy import stats
from scipy.stats import norm
from scipy.optimize import curve_fit

# ファイル読み込みのためのやつ
import re

# webに飛ぶやつ
import webbrowser

#windowの作成
root=tk.Tk()
#windowのタイトル変更
root.title(f"TAS simulation ver: {__version__}")
#windowのサイズ指定
root.geometry("550x840")#550*840

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

# GUIの配分を決める。
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=2)
root.rowconfigure(1, weight=2)
root.rowconfigure(2, weight=4)
root.rowconfigure(3, weight=2)
root.rowconfigure(4, weight=5)
root.rowconfigure(5, weight=3)

# ファイル選択のフレームの作成と設置
frame1 = ttk.Labelframe(root,text= "lattice infromation")
frame1.grid(row=0,column=0,sticky="NSEW")

frame1.columnconfigure(0, weight=1)
frame1.columnconfigure(1, weight=1)
frame1.columnconfigure(2, weight=1)
frame1.columnconfigure(3, weight=1)
frame1.columnconfigure(4, weight=1)
frame1.columnconfigure(5, weight=1)
frame1.rowconfigure(0, weight=1)
frame1.rowconfigure(1, weight=1)

from RLcalc import RL_calc  #
from UBcalc import UB_calc  #
from anglecalc import angle_calc    #

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
la.insert(0,'5.026307')

lc2 = tk.Label(frame1,text='b (Å)')
lc2.grid(row=0, column=1,sticky="NSEW")
lb = ttk.Entry(frame1)
lb.grid(row=1, column=1,sticky="NSEW")
lb.insert(0,'5.026307')

lc3 = tk.Label(frame1,text='c (Å)')
lc3.grid(row=0, column=2,sticky="NSEW")
lc = ttk.Entry(frame1)
lc.grid(row=1, column=2,sticky="NSEW")
lc.insert(0,'13.784500')

lc4 = tk.Label(frame1,text='α (deg)')
lc4.grid(row=0, column=3,sticky="NSEW")
lc_alpha = ttk.Entry(frame1)
lc_alpha.grid(row=1, column=3,sticky="NSEW")
lc_alpha.insert(0,'90')

lc5 = tk.Label(frame1,text='β (deg)')
lc5.grid(row=0, column=4,sticky="NSEW")
lc_beta = ttk.Entry(frame1)
lc_beta.grid(row=1, column=4,sticky="NSEW")
lc_beta.insert(0,'90')

lc6 = tk.Label(frame1,text='γ (deg)')
lc6.grid(row=0, column=5,sticky="NSEW")
lc_gamma = ttk.Entry(frame1)
lc_gamma.grid(row=1, column=5,sticky="NSEW")
lc_gamma.insert(0,'120')

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
sv1_h.insert(0,'1')

sv1 = tk.Label(frame2a,text='k')
sv1.grid(row=0, column=1,sticky="NSEW")
sv1_k = ttk.Entry(frame2a)
sv1_k.grid(row=1, column=1,sticky="NSEW")
sv1_k.insert(0,'0')

sv1 = tk.Label(frame2a,text='l')
sv1.grid(row=0, column=2,sticky="NSEW")
sv1_l = ttk.Entry(frame2a)
sv1_l.grid(row=1, column=2,sticky="NSEW")
sv1_l.insert(0,'0')

sv2 = tk.Label(frame2b,text='h')
sv2.grid(row=0, column=0,sticky="NSEW")
sv2_h = ttk.Entry(frame2b)
sv2_h.grid(row=1, column=0,sticky="NSEW")
sv2_h.insert(0,'0')

sv2 = tk.Label(frame2b,text='k')
sv2.grid(row=0, column=1,sticky="NSEW")
sv2_k = ttk.Entry(frame2b)
sv2_k.grid(row=1, column=1,sticky="NSEW")
sv2_k.insert(0,'1')

sv2 = tk.Label(frame2b,text='l')
sv2.grid(row=0, column=2,sticky="NSEW")
sv2_l = ttk.Entry(frame2b)
sv2_l.grid(row=1, column=2,sticky="NSEW")
sv2_l.insert(0,'0')

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
def calculate_all():
    # on_Rlcalc の結果を取得
    RLtable = on_Rlcalc()
    #print("RLtable:", RLtable)
    
    # on_UBcalc の結果を取得
    UBtable = on_UBcalc()
    #print("UBtable:", UBtable)
    
    # U matrixを表示
    txt_u_11.config(state="normal")  # 一時的に編集可能に
    txt_u_11.delete(0, tk.END)
    txt_u_11.insert(0, round(UBtable['UB'][0,0],4))
    txt_u_11.config(state="readonly") # 編集不可に設定
    txt_u_12.config(state="normal")  # 一時的に編集可能に
    txt_u_12.delete(0, tk.END)
    txt_u_12.insert(0, round(UBtable['UB'][0,1],4))
    txt_u_12.config(state="readonly") # 編集不可に設定
    txt_u_13.config(state="normal")  # 一時的に編集可能に
    txt_u_13.delete(0, tk.END)
    txt_u_13.insert(0, round(UBtable['UB'][0,2],4))
    txt_u_13.config(state="readonly") # 編集不可に設定
    txt_u_21.config(state="normal")  # 一時的に編集可能に
    txt_u_21.delete(0, tk.END)
    txt_u_21.insert(0, round(UBtable['UB'][1,0],4))
    txt_u_21.config(state="readonly") # 編集不可に設定
    txt_u_22.config(state="normal")  # 一時的に編集可能に
    txt_u_22.delete(0, tk.END)
    txt_u_22.insert(0, round(UBtable['UB'][1,1],4))
    txt_u_22.config(state="readonly") # 編集不可に設定
    txt_u_23.config(state="normal")  # 一時的に編集可能に
    txt_u_23.delete(0, tk.END)
    txt_u_23.insert(0, round(UBtable['UB'][1,2],4))
    txt_u_23.config(state="readonly") # 編集不可に設定
    txt_u_31.config(state="normal")  # 一時的に編集可能に
    txt_u_31.delete(0, tk.END)
    txt_u_31.insert(0, round(UBtable['UB'][2,0],4))
    txt_u_31.config(state="readonly") # 編集不可に設定
    txt_u_32.config(state="normal")  # 一時的に編集可能に
    txt_u_32.delete(0, tk.END)
    txt_u_32.insert(0, round(UBtable['UB'][2,1],4))
    txt_u_32.config(state="readonly") # 編集不可に設定
    txt_u_33.config(state="normal")  # 一時的に編集可能に
    txt_u_33.delete(0, tk.END)
    txt_u_33.insert(0, round(UBtable['UB'][2,2],4))
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

#ボタン1つで両方の計算を実行
UBcalculate_button = tk.Button(frame3, text="UB calculation", command=calculate_all)
UBcalculate_button.grid(row=1, column=1,sticky="NSEW")

# ブラッグピーク位置を入力
# ファイル選択のフレームの作成と設置
frame4 = ttk.Labelframe(root,text= "bragg peak position")
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
eief.set(1)

# ラジオボタンを作成し、commandにupdate_labelを設定
rdo_eief0 = tk.Radiobutton(frame4, value=0, variable=eief, text='Ei fix', command=update_label, width=15)
rdo_eief0.grid(row=0, column=0, sticky="NSEW")

rdo_eief1 = tk.Radiobutton(frame4, value=1, variable=eief, text='Ef fix', command=update_label, width=15)
rdo_eief1.grid(row=1, column=0, sticky="NSEW")

bpl1 = tk.Label(frame4,text='Ef')
bpl1.grid(row=0, column=1,sticky="NSEW")
Energy = ttk.Entry(frame4)
Energy.grid(row=1, column=1,sticky="NSEW")
Energy.insert(0,'4.8')

bpl2 = tk.Label(frame4,text='h')
bpl2.grid(row=0, column=2,sticky="NSEW")
bp_h = ttk.Entry(frame4)
bp_h.grid(row=1, column=2,sticky="NSEW")
bp_h.insert(0,'1.5')

bpl3 = tk.Label(frame4,text='k')
bpl3.grid(row=0, column=3,sticky="NSEW")
bp_k = ttk.Entry(frame4)
bp_k.grid(row=1, column=3,sticky="NSEW")
bp_k.insert(0,'0')

bpl4 = tk.Label(frame4,text='l')
bpl4.grid(row=0, column=4,sticky="NSEW")
bp_l = ttk.Entry(frame4)
bp_l.grid(row=1, column=4,sticky="NSEW")
bp_l.insert(0,'0')

bpl5 = tk.Label(frame4,text='C2')
bpl5.grid(row=0, column=5,sticky="NSEW")
bp_c2 = ttk.Entry(frame4)
bp_c2.grid(row=1, column=5,sticky="NSEW")
bp_c2.insert(0,'49.75')

bpl6 = tk.Label(frame4,text='μ')
bpl6.grid(row=0, column=6,sticky="NSEW")
bp_mu = ttk.Entry(frame4)
bp_mu.grid(row=1, column=6,sticky="NSEW")
bp_mu.insert(0,'0')

bpl7 = tk.Label(frame4,text='ν')
bpl7.grid(row=0, column=7,sticky="NSEW")
bp_nu = ttk.Entry(frame4)
bp_nu.grid(row=1, column=7,sticky="NSEW")
bp_nu.insert(0,'0')

# select scan type
# ファイル選択のフレームの作成と設置
frame6 = ttk.Labelframe(root,text= "select scan type")
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
tab_001 = tk.Frame(notebook00)# Angle calculation
tab_002 = tk.Frame(notebook00)# constant Q scan
tab_003 = tk.Frame(notebook00)# constant E scan
# notebookにタブを追加
notebook00.add(tab_001, text="single Q-E point")
notebook00.add(tab_002, text="constant Q scan")
notebook00.add(tab_003, text="constant E scan")

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
    
    # 内積の計算
    dat_ax1_ax2 = np.dot(axis1, axis2)
    
    # コサインの逆関数を用いて角度を計算 (ラジアン)
    angle_ax1_ax2 = np.degrees(np.arccos(dat_ax1_ax2 / (norm_sv1 * norm_sv2)))
    
    U=UBtable['U']
    B=UBtable['B']
    UB=UBtable['UB']
    
    # Bragg peak positionの取得
    bpe = float(Energy.get())
    bpc2 = float(bp_c2.get())
    bpmu = float(bp_mu.get())
    bpnu = float(bp_nu.get())
    bp = np.array([float(bp_h.get()), float(bp_k.get()), float(bp_l.get())])
    
    # calculation pointの取得
    cphw = float(acbe.get())
    cp = np.array([float(acbh.get()), float(acbk.get()), float(acbl.get())])
    
    # Ei or Ef fixの判定
    fixe=float(eief.get())
    
    # UBtableを計算
    angletable = angle_calc(
        astar,bstar,cstar,UB,bpe,bpc2,bpmu,bpnu,bp,cphw,cp,fixe,angle_ax1_ax2
    )
    
    return angletable

def calculate_angle():
    # on_anglecalc の結果を取得
    angletable = on_anglecalc()
    #print("RLtable:", RLtable)
    
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
Angle_calculate_button = tk.Button(tab_001a, text="Angle calculation", command=calculate_angle)
Angle_calculate_button.grid(row=1, column=4,sticky="NSEW")

tab_001b = ttk.Labelframe(tab_001,text= "calculation results")
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

acl2 = tk.Label(tab_001b,text='A1')
acl2.grid(row=0, column=1,sticky="NSEW")
acb2 = ttk.Entry(tab_001b)
acb2.grid(row=1, column=1,sticky="NSEW")
acb2.insert(0,'0')

acl3 = tk.Label(tab_001b,text='C2')
acl3.grid(row=0, column=2,sticky="NSEW")
acb3 = ttk.Entry(tab_001b)
acb3.grid(row=1, column=2,sticky="NSEW")
acb3.insert(0,'0')

acl4 = tk.Label(tab_001b,text='A2')
acl4.grid(row=0, column=3,sticky="NSEW")
acb4 = ttk.Entry(tab_001b)
acb4.grid(row=1, column=3,sticky="NSEW")
acb4.insert(0,'0')

acl5 = tk.Label(tab_001b,text='C3')
acl5.grid(row=0, column=4,sticky="NSEW")
acb5 = ttk.Entry(tab_001b)
acb5.grid(row=1, column=4,sticky="NSEW")
acb5.insert(0,'0')

acl6 = tk.Label(tab_001b,text='A3')
acl6.grid(row=0, column=5,sticky="NSEW")
acb6 = ttk.Entry(tab_001b)
acb6.grid(row=1, column=5,sticky="NSEW")
acb6.insert(0,'0')

acl7 = tk.Label(tab_001b,text='μ')
acl7.grid(row=0, column=6,sticky="NSEW")
acb7 = ttk.Entry(tab_001b)
acb7.grid(row=1, column=6,sticky="NSEW")
acb7.insert(0,'0')

acl8 = tk.Label(tab_001b,text='ν')
acl8.grid(row=0, column=7,sticky="NSEW")
acb8 = ttk.Entry(tab_001b)
acb8.grid(row=1, column=7,sticky="NSEW")
acb8.insert(0,'0')

acl9 = tk.Label(tab_001b,text='warning : ')
acl9.grid(row=2, column=0,sticky="NSEW")

# hardware limit
# ファイル選択のフレームの作成と設置
frame5 = ttk.Labelframe(root,text= "hardware limit")
frame5.grid(row=5,column=0,sticky="NSEW")

frame5.columnconfigure(0, weight=1)
frame5.columnconfigure(1, weight=1)
frame5.columnconfigure(2, weight=1)
frame5.columnconfigure(3, weight=1)
frame5.columnconfigure(4, weight=1)
frame5.columnconfigure(5, weight=1)
frame5.columnconfigure(6, weight=1)
frame5.columnconfigure(7, weight=1)
frame5.columnconfigure(8, weight=1)
frame5.rowconfigure(0, weight=1)
frame5.rowconfigure(1, weight=1)
frame5.rowconfigure(2, weight=1)

hwll0 = tk.Label(frame5,text='from',width=15)
hwll0.grid(row=1, column=0,sticky="NSEW")
hwll1 = tk.Label(frame5,text='to',width=15)
hwll1.grid(row=2, column=0,sticky="NSEW")

hwll2 = tk.Label(frame5,text='C1')
hwll2.grid(row=0, column=1,sticky="NSEW")
hwl2f = ttk.Entry(frame5)
hwl2f.grid(row=1, column=1,sticky="NSEW")
hwl2f.insert(0,'19.9305')
hwl2t = ttk.Entry(frame5)
hwl2t.grid(row=2, column=1,sticky="NSEW")
hwl2t.insert(0,'58.482')

hwll3 = tk.Label(frame5,text='A1')
hwll3.grid(row=0, column=2,sticky="NSEW")
hwl3f = ttk.Entry(frame5)
hwl3f.grid(row=1, column=2,sticky="NSEW")
hwl3f.insert(0,'39.861')
hwl3t = ttk.Entry(frame5)
hwl3t.grid(row=2, column=2,sticky="NSEW")
hwl3t.insert(0,'116.964')

hwll4 = tk.Label(frame5,text='C2')
hwll4.grid(row=0, column=3,sticky="NSEW")
hwl4f = ttk.Entry(frame5)
hwl4f.grid(row=1, column=3,sticky="NSEW")
hwl4f.insert(0,'-180')
hwl4t = ttk.Entry(frame5)
hwl4t.grid(row=2, column=3,sticky="NSEW")
hwl4t.insert(0,'180')

hwll5 = tk.Label(frame5,text='A2')
hwll5.grid(row=0, column=4,sticky="NSEW")
hwl5f = ttk.Entry(frame5)
hwl5f.grid(row=1, column=4,sticky="NSEW")
hwl5f.insert(0,'6')
hwl5t = ttk.Entry(frame5)
hwl5t.grid(row=2, column=4,sticky="NSEW")
hwl5t.insert(0,'120')

hwll6 = tk.Label(frame5,text='C3')
hwll6.grid(row=0, column=5,sticky="NSEW")
hwl6f = ttk.Entry(frame5)
hwl6f.grid(row=1, column=5,sticky="NSEW")
hwl6f.insert(0,'19.9305')
hwl6t = ttk.Entry(frame5)
hwl6t.grid(row=2, column=5,sticky="NSEW")
hwl6t.insert(0,'58.482')

hwll7 = tk.Label(frame5,text='A3')
hwll7.grid(row=0, column=6,sticky="NSEW")
hwl7f = ttk.Entry(frame5)
hwl7f.grid(row=1, column=6,sticky="NSEW")
hwl7f.insert(0,'39.861')
hwl7t = ttk.Entry(frame5)
hwl7t.grid(row=2, column=6,sticky="NSEW")
hwl7t.insert(0,'116.964')

hwll8 = tk.Label(frame5,text='μ')
hwll8.grid(row=0, column=7,sticky="NSEW")
hwl8f = ttk.Entry(frame5)
hwl8f.grid(row=1, column=7,sticky="NSEW")
hwl8f.insert(0,'-5')
hwl8t = ttk.Entry(frame5)
hwl8t.grid(row=2, column=7,sticky="NSEW")
hwl8t.insert(0,'5')

hwll9 = tk.Label(frame5,text='ν')
hwll9.grid(row=0, column=8,sticky="NSEW")
hwl9f = ttk.Entry(frame5)
hwl9f.grid(row=1, column=8,sticky="NSEW")
hwl9f.insert(0,'-5')
hwl9t = ttk.Entry(frame5)
hwl9t.grid(row=2, column=8,sticky="NSEW")
hwl9t.insert(0,'5')

#window状態の維持
root.mainloop()

#############
# pyinstaller tips
# pyinstaller -F --noconsole display.py
"""
コマンド
pyinstaller ASYURA.py --noconsole
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