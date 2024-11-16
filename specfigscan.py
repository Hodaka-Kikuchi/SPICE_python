import configparser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
from matplotlib.widgets import Slider
import os
import sys

def plot_spectrometer(A_sets,QE_sets, initial_index=0):
    """
    三軸分光器のXY平面図を描く
    :param A_sets: A1, A2, A3 の角度セットのリスト
    :param initial_index: 初期表示する角度セットのインデックス
    """
    # 距離と半径の設定
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
    
    # 設定を変数に代入
    monochromator_radius = int(config['settings']['monochromator_radius'])
    monochromator_to_sample = int(config['settings']['monochromator_to_sample'])
    sample_goniometer_radius = int(config['settings']['sample_goniometer_radius'])
    sample_to_analyzer = int(config['settings']['sample_to_analyzer'])
    analyzer_radius = int(config['settings']['analyzer_radius'])
    analyzer_to_detector = int(config['settings']['analyzer_to_detector'])
    detector_radius = int(config['settings']['detector_radius'])
    floor_length = int(config['settings']['floor_length'])
    floor_width = int(config['settings']['floor_width'])
    floor_position_x = int(config['settings']['floor_position_x'])
    floor_position_y = int(config['settings']['floor_position_y'])

    # 距離のマックス
    max_L = monochromator_to_sample + sample_to_analyzer + analyzer_to_detector

    # モノクロメータの位置 (原点)
    x_mono = 0
    y_mono = 0

    # プロット設定
    fig, ax = plt.subplots(figsize=(8, 5))  # 図のサイズを縮小
    plt.subplots_adjust(left=0.01, bottom=0.25)  # 左の余白を削る

    # ソース
    ax.scatter(-max_L, 0, color='black', s=1, label='Source')

    # スライダー設定
    ax_slider = plt.axes([0.25, 0.10, 0.65, 0.03], facecolor='lightgoldenrodyellow')
    slider = Slider(ax_slider, 'scan number', 1, len(A_sets) , valinit=initial_index+1, valstep=1)

    # スキャン条件表示用テキストを初期化
    ax.text(0.4, 1.05, f'ℏω: {QE_sets[initial_index][0]} meV, h: {QE_sets[initial_index][1]}, k: {QE_sets[initial_index][2]}, l: {QE_sets[initial_index][3]}', 
                                   horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

    def update(val):
        # スライダーの値に基づいて角度セットを選択
        index = int(val)-1
        A1, A2, A3 = A_sets[index]

        # 位置計算
        x_sample = monochromator_to_sample * np.cos(np.radians(A1))
        y_sample = monochromator_to_sample * np.sin(np.radians(A1))
        x_analyzer = x_sample + sample_to_analyzer * np.cos(np.radians(A1 + A2))
        y_analyzer = y_sample + sample_to_analyzer * np.sin(np.radians(A1 + A2))
        x_detector = x_analyzer + analyzer_to_detector * np.cos(np.radians(A1 + A2 + A3))
        y_detector = y_analyzer + analyzer_to_detector * np.sin(np.radians(A1 + A2 + A3))

        # プロットを再描画
        ax.clear()
        ax.scatter(-max_L, 0, color='black', s=1, label='Source')
        ax.add_patch(patches.Circle((x_mono, y_mono), monochromator_radius, linewidth=2, edgecolor='none', facecolor='blue', linestyle='--'))
        ax.scatter(x_mono, y_mono, color='blue', s=monochromator_radius, label='Monochromator')

        ax.add_patch(patches.Circle((x_sample, y_sample), sample_goniometer_radius, linewidth=2, edgecolor='none', facecolor='green', linestyle='--'))
        ax.scatter(x_sample, y_sample, color='green', s=sample_goniometer_radius, label='Sample Gonio')

        ax.add_patch(patches.Circle((x_analyzer, y_analyzer), analyzer_radius, linewidth=2, edgecolor='none', facecolor='red', linestyle='--'))
        ax.scatter(x_analyzer, y_analyzer, color='red', s=analyzer_radius, label='Analyzer')

        ax.add_patch(patches.Circle((x_detector, y_detector), detector_radius, linewidth=2, edgecolor='none', facecolor='purple', linestyle='--'))
        ax.scatter(x_detector, y_detector, color='purple', s=detector_radius, label='Detector')

        # 各コンポーネント間の線を引く
        ax.plot([-max_L, x_mono], [0, y_mono], color='black', linestyle='--')
        ax.plot([x_mono, x_sample], [y_mono, y_sample], color='blue', linestyle='--')
        ax.plot([x_sample, x_analyzer], [y_sample, y_analyzer], color='green', linestyle='--')
        ax.plot([x_analyzer, x_detector], [y_analyzer, y_detector], color='red', linestyle='--')

        # floorの外形を書く
        ax.add_patch(Rectangle((floor_position_x, floor_position_y), floor_length, floor_width, linewidth=1, edgecolor='black', facecolor='none', linestyle='--'))

        # 軸設定
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim(-max_L, max_L)
        ax.set_ylim(-800, max_L * 1.25)

        ax.set_xlabel('x [m]')  # 横軸ラベル
        ax.set_ylabel('y [m]')  # 縦軸ラベル
        
        ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
        
        # スキャン条件の更新
        ax.text(0.4, 1.05, f'ℏω: {QE_sets[index][0]} meV, h: {QE_sets[index][1]}, k: {QE_sets[index][2]}, l: {QE_sets[index][3]}', 
                                   horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        plt.draw()

    slider.on_changed(update)

    # 初期の描画を行う
    update(slider.val)

    def on_key(event):
        """
        左右矢印キーでスライダーを動かす
        """
        if event.key == 'right':
            slider.set_val(min(slider.val + 1, len(A_sets)))
        elif event.key == 'left':
            slider.set_val(max(slider.val - 1, 1))

    # キーイベントを設定
    fig.canvas.mpl_connect('key_press_event', on_key)

    plt.show()
