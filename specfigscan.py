import configparser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
from matplotlib.widgets import Slider
import os
import sys

def plot_spectrometer(inst_param,sense,A_sets,QE_sets, initial_index=0):
    """
    三軸分光器のXY平面図を描く
    :param A_sets: A1, A2, A3 の角度セットのリスト
    :param initial_index: 初期表示する角度セットのインデックス
    """
    '''
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
    '''
    # 設定を変数に代入
    # systemセクションのview設定を読み込む
    mono_radius = inst_param["mono_radius"]
    mono_to_sample = inst_param["mono_to_sample"]
    sample_stage_radius = inst_param["sample_stage_radius"]
    sample_to_ana = inst_param["sample_to_ana"]
    ana_radius = inst_param["ana_radius"]
    ana_to_det = inst_param["ana_to_det"]
    det_radius = inst_param["det_radius"]
    floor_length = inst_param["floor_length"]
    floor_width = inst_param["floor_width"]
    floor_position_x = inst_param["floor_position_x"]
    floor_position_y = inst_param["floor_position_y"]

    # 距離のマックス
    max_L = mono_to_sample + sample_to_ana + ana_to_det

    # モノクロメータの位置 (原点)
    x_mono = 0
    y_mono = 0

    # プロット設定
    fig, ax = plt.subplots(figsize=(8, 5))  # 図のサイズを縮小
    plt.subplots_adjust(left=0.01, bottom=0.25)  # 左の余白を削る

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
        x_sample = mono_to_sample * np.cos(np.radians(A1))
        y_sample = mono_to_sample * np.sin(np.radians(A1))
        x_ana = x_sample + sample_to_ana * np.cos(np.radians(A1 + A2))
        y_ana = y_sample + sample_to_ana * np.sin(np.radians(A1 + A2))
        x_det = x_ana + ana_to_det * np.cos(np.radians(A1 + A2 + A3))
        y_det = y_ana + ana_to_det * np.sin(np.radians(A1 + A2 + A3))

        # プロットを再描画
        ax.clear()
        ax.scatter(-max_L, 0, color='black', s=10, label='Source')
        ax.add_patch(patches.Circle((x_mono, y_mono), mono_radius, linewidth=2, edgecolor='none', facecolor='blue', linestyle='--'))
        ax.scatter(x_mono, y_mono, color='blue', s=10, label='mono')

        ax.add_patch(patches.Circle((x_sample, y_sample), sample_stage_radius, linewidth=2, edgecolor='none', facecolor='green', linestyle='--'))
        ax.scatter(x_sample, y_sample, color='green', s=10, label='Sample Gonio')

        ax.add_patch(patches.Circle((x_ana, y_ana), ana_radius, linewidth=2, edgecolor='none', facecolor='red', linestyle='--'))
        ax.scatter(x_ana, y_ana, color='red', s=10, label='ana')

        ax.add_patch(patches.Circle((x_det, y_det), det_radius, linewidth=2, edgecolor='none', facecolor='purple', linestyle='--'))
        ax.scatter(x_det, y_det, color='purple', s=10, label='det')

        # 各コンポーネント間の線を引く
        ax.plot([-max_L, x_mono], [0, y_mono], color='black', linestyle='--')
        ax.plot([x_mono, x_sample], [y_mono, y_sample], color='blue', linestyle='--')
        ax.plot([x_sample, x_ana], [y_sample, y_ana], color='green', linestyle='--')
        ax.plot([x_ana, x_det], [y_ana, y_det], color='red', linestyle='--')

        # floorの外形を書く
        ax.add_patch(Rectangle((floor_position_x, floor_position_y), floor_length, floor_width, linewidth=1, edgecolor='black', facecolor='none', linestyle='--'))

        # 軸設定
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim(-floor_length/2 * 1.1, floor_length/2 * 1.1)
        
        # system設定に基づいてy軸の設定を変更
        if sense == '-+-':
            ax.set_ylim((floor_width+mono_radius) * 1.1, -mono_radius * 1.1)  # 上下反転
        elif sense == '+-+':
            ax.set_ylim(-mono_radius * 1.1, (floor_width+mono_radius) * 1.1)  # 通常の範囲

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

    #plt.show()
