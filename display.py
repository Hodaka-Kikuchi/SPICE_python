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
root.rowconfigure(0, weight=1)
root.rowconfigure(1, weight=1)
root.rowconfigure(2, weight=1)
root.rowconfigure(3, weight=1)
root.rowconfigure(4, weight=1)
root.rowconfigure(5, weight=1)
root.rowconfigure(6, weight=1)

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
la = ttk.Entry(frame1,width=5)
la.grid(row=1, column=0,sticky="NSEW")
la.insert(0,'7.1963')

lc2 = tk.Label(frame1,text='b (Å)')
lc2.grid(row=0, column=1,sticky="NSEW")
lb = ttk.Entry(frame1,width=5)
lb.grid(row=1, column=1,sticky="NSEW")
lb.insert(0,'7.1963')

lc3 = tk.Label(frame1,text='c (Å)')
lc3.grid(row=0, column=2,sticky="NSEW")
lc = ttk.Entry(frame1,width=5)
lc.grid(row=1, column=2,sticky="NSEW")
lc.insert(0,'6')

lc4 = tk.Label(frame1,text='α (deg)')
lc4.grid(row=0, column=3,sticky="NSEW")
lc_alpha = ttk.Entry(frame1,width=5)
lc_alpha.grid(row=1, column=3,sticky="NSEW")
lc_alpha.insert(0,'90')

lc5 = tk.Label(frame1,text='β (deg)')
lc5.grid(row=0, column=4,sticky="NSEW")
lc_beta = ttk.Entry(frame1,width=5)
lc_beta.grid(row=1, column=4,sticky="NSEW")
lc_beta.insert(0,'90')

lc6 = tk.Label(frame1,text='γ (deg)')
lc6.grid(row=0, column=5,sticky="NSEW")
lc_gamma = ttk.Entry(frame1,width=5)
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
sv1_h = ttk.Entry(frame2a,width=5)
sv1_h.grid(row=1, column=0,sticky="NSEW")
sv1_h.insert(0,'1')

sv1 = tk.Label(frame2a,text='k')
sv1.grid(row=0, column=1,sticky="NSEW")
sv1_k = ttk.Entry(frame2a,width=5)
sv1_k.grid(row=1, column=1,sticky="NSEW")
sv1_k.insert(0,'0')

sv1 = tk.Label(frame2a,text='l')
sv1.grid(row=0, column=2,sticky="NSEW")
sv1_l = ttk.Entry(frame2a,width=5)
sv1_l.grid(row=1, column=2,sticky="NSEW")
sv1_l.insert(0,'0')

sv2 = tk.Label(frame2b,text='h')
sv2.grid(row=0, column=0,sticky="NSEW")
sv2_h = ttk.Entry(frame2b,width=5)
sv2_h.grid(row=1, column=0,sticky="NSEW")
sv2_h.insert(0,'0')

sv2 = tk.Label(frame2b,text='k')
sv2.grid(row=0, column=1,sticky="NSEW")
sv2_k = ttk.Entry(frame2b,width=5)
sv2_k.grid(row=1, column=1,sticky="NSEW")
sv2_k.insert(0,'0')

sv2 = tk.Label(frame2b,text='l')
sv2.grid(row=0, column=2,sticky="NSEW")
sv2_l = ttk.Entry(frame2b,width=5)
sv2_l.grid(row=1, column=2,sticky="NSEW")
sv2_l.insert(0,'3')

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
txt_u_11 = ttk.Entry(frame3a,width=5,state="readonly")
txt_u_11.insert(0,'0')
txt_u_11.grid(row=0, column=0,sticky="NSEW")
txt_u_12 = ttk.Entry(frame3a,width=5,state="readonly")
txt_u_12.insert(0,'0')
txt_u_12.grid(row=0, column=1,sticky="NSEW")
txt_u_13 = ttk.Entry(frame3a,width=5,state="readonly")
txt_u_13.insert(0,'0')
txt_u_13.grid(row=0, column=2,sticky="NSEW")
txt_u_21 = ttk.Entry(frame3a,width=5,state="readonly")
txt_u_21.insert(0,'0')
txt_u_21.grid(row=1, column=0,sticky="NSEW")
txt_u_22 = ttk.Entry(frame3a,width=5,state="readonly")
txt_u_22.insert(0,'0')
txt_u_22.grid(row=1, column=1,sticky="NSEW")
txt_u_23 = ttk.Entry(frame3a,width=5,state="readonly")
txt_u_23.insert(0,'0')
txt_u_23.grid(row=1, column=2,sticky="NSEW")
txt_u_31 = ttk.Entry(frame3a,width=5,state="readonly")
txt_u_31.insert(0,'0')
txt_u_31.grid(row=2, column=0,sticky="NSEW")
txt_u_32 = ttk.Entry(frame3a,width=5,state="readonly")
txt_u_32.insert(0,'0')
txt_u_32.grid(row=2, column=1,sticky="NSEW")
txt_u_33 = ttk.Entry(frame3a,width=5,state="readonly")
txt_u_33.insert(0,'0')
txt_u_33.grid(row=2, column=2,sticky="NSEW")

# Bmatrixの表示
txt_b_11 = ttk.Entry(frame3b,width=5,state="readonly")
txt_b_11.insert(0,'0')
txt_b_11.grid(row=0, column=0,sticky="NSEW")
txt_b_12 = ttk.Entry(frame3b,width=5,state="readonly")
txt_b_12.insert(0,'0')
txt_b_12.grid(row=0, column=1,sticky="NSEW")
txt_b_13 = ttk.Entry(frame3b,width=5,state="readonly")
txt_b_13.insert(0,'0')
txt_b_13.grid(row=0, column=2,sticky="NSEW")
txt_b_21 = ttk.Entry(frame3b,width=5,state="readonly")
txt_b_21.insert(0,'0')
txt_b_21.grid(row=1, column=0,sticky="NSEW")
txt_b_22 = ttk.Entry(frame3b,width=5,state="readonly")
txt_b_22.insert(0,'0')
txt_b_22.grid(row=1, column=1,sticky="NSEW")
txt_b_23 = ttk.Entry(frame3b,width=5,state="readonly")
txt_b_23.insert(0,'0')
txt_b_23.grid(row=1, column=2,sticky="NSEW")
txt_b_31 = ttk.Entry(frame3b,width=5,state="readonly")
txt_b_31.insert(0,'0')
txt_b_31.grid(row=2, column=0,sticky="NSEW")
txt_b_32 = ttk.Entry(frame3b,width=5,state="readonly")
txt_b_32.insert(0,'0')
txt_b_32.grid(row=2, column=1,sticky="NSEW")
txt_b_33 = ttk.Entry(frame3b,width=5,state="readonly")
txt_b_33.insert(0,'0')
txt_b_33.grid(row=2, column=2,sticky="NSEW")

# UBmatrixの表示
txt_ub_11 = ttk.Entry(frame3c,width=5,state="readonly")
txt_ub_11.insert(0,'0')
txt_ub_11.grid(row=0, column=0,sticky="NSEW")
txt_ub_12 = ttk.Entry(frame3c,width=5,state="readonly")
txt_ub_12.insert(0,'0')
txt_ub_12.grid(row=0, column=1,sticky="NSEW")
txt_ub_13 = ttk.Entry(frame3c,width=5,state="readonly")
txt_ub_13.insert(0,'0')
txt_ub_13.grid(row=0, column=2,sticky="NSEW")
txt_ub_21 = ttk.Entry(frame3c,width=5,state="readonly")
txt_ub_21.insert(0,'0')
txt_ub_21.grid(row=1, column=0,sticky="NSEW")
txt_ub_22 = ttk.Entry(frame3c,width=5,state="readonly")
txt_ub_22.insert(0,'0')
txt_ub_22.grid(row=1, column=1,sticky="NSEW")
txt_ub_23 = ttk.Entry(frame3c,width=5,state="readonly")
txt_ub_23.insert(0,'0')
txt_ub_23.grid(row=1, column=2,sticky="NSEW")
txt_ub_31 = ttk.Entry(frame3c,width=5,state="readonly")
txt_ub_31.insert(0,'0')
txt_ub_31.grid(row=2, column=0,sticky="NSEW")
txt_ub_32 = ttk.Entry(frame3c,width=5,state="readonly")
txt_ub_32.insert(0,'0')
txt_ub_32.grid(row=2, column=1,sticky="NSEW")
txt_ub_33 = ttk.Entry(frame3c,width=5,state="readonly")
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
calculate_button = tk.Button(frame3, text="UB calculation", command=calculate_all)
calculate_button.grid(row=1, column=1,sticky="NSEW")

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
frame4.rowconfigure(0, weight=1)
frame4.rowconfigure(1, weight=1)

EfEi = tk.Label(frame4,text='Ef (meV)')
EfEi.grid(row=0, column=0,sticky="NSEW")
Energy = ttk.Entry(frame4,width=5)
Energy.grid(row=1, column=0,sticky="NSEW")
Energy.insert(0,'5')

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