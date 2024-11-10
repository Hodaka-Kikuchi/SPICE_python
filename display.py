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
frame1.rowconfigure(2, weight=1)
frame1.rowconfigure(3, weight=1)

from UBcalc import UB_calc  #
#from RLcalc import calc_RL  #

def on_UBcalc():
    a = float(la.get())
    b = float(lb.get())
    c = float(lc.get())
    alpha = float(lc_alpha.get())
    beta = float(lc_beta.get())
    gamma = float(lc_gamma.get())

    result = UB_calc(a,b,c,alpha,beta,gamma)
    #return astar,bstar,cstar,alpha_star,beta_star,gamma_star
    print(result)  
    
calculate_button = tk.Button(frame1, text="計算", command=on_UBcalc)
calculate_button.grid(row=4, column=0, columnspan=6)

# 格子定数を入力する欄
lc1 = tk.Label(frame1,text='a (Å)')
lc1.grid(row=0, column=0,sticky="NSEW")
la = ttk.Entry(frame1,width=5)
la.grid(row=1, column=0,sticky="NSEW")

lc2 = tk.Label(frame1,text='b (Å)')
lc2.grid(row=0, column=1,sticky="NSEW")
lb = ttk.Entry(frame1,width=5)
lb.grid(row=1, column=1,sticky="NSEW")

lc3 = tk.Label(frame1,text='c (Å)')
lc3.grid(row=0, column=2,sticky="NSEW")
lc = ttk.Entry(frame1,width=5)
lc.grid(row=1, column=2,sticky="NSEW")

lc4 = tk.Label(frame1,text='α (deg)')
lc4.grid(row=0, column=3,sticky="NSEW")
lc_alpha = ttk.Entry(frame1,width=5)
lc_alpha.grid(row=1, column=3,sticky="NSEW")

lc5 = tk.Label(frame1,text='β (deg)')
lc5.grid(row=0, column=4,sticky="NSEW")
lc_beta = ttk.Entry(frame1,width=5)
lc_beta.grid(row=1, column=4,sticky="NSEW")

lc6 = tk.Label(frame1,text='γ (deg)')
lc6.grid(row=0, column=5,sticky="NSEW")
lc_gamma = ttk.Entry(frame1,width=5)
lc_gamma.grid(row=1, column=5,sticky="NSEW")

# 散乱面を入力する欄
sp1 = tk.Label(frame1,text='h')
sp1.grid(row=2, column=0,sticky="NSEW")
sp1_h = ttk.Entry(frame1,width=5)
sp1_h.grid(row=3, column=0,sticky="NSEW")

sp1 = tk.Label(frame1,text='k')
sp1.grid(row=2, column=1,sticky="NSEW")
sp1_k = ttk.Entry(frame1,width=5)
sp1_k.grid(row=3, column=1,sticky="NSEW")

sp1 = tk.Label(frame1,text='l')
sp1.grid(row=2, column=2,sticky="NSEW")
sp1_l = ttk.Entry(frame1,width=5)
sp1_l.grid(row=3, column=2,sticky="NSEW")

sp2 = tk.Label(frame1,text='h')
sp2.grid(row=2, column=3,sticky="NSEW")
sp2_h = ttk.Entry(frame1,width=5)
sp2_h.grid(row=3, column=3,sticky="NSEW")

sp2 = tk.Label(frame1,text='k')
sp2.grid(row=2, column=4,sticky="NSEW")
sp2_k = ttk.Entry(frame1,width=5)
sp2_k.grid(row=3, column=4,sticky="NSEW")

sp2 = tk.Label(frame1,text='l')
sp2.grid(row=2, column=5,sticky="NSEW")
sp2_l = ttk.Entry(frame1,width=5)
sp2_l.grid(row=3, column=5,sticky="NSEW")

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