import configparser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
from matplotlib.widgets import Slider
from PIL import Image  # GIF 保存のために必要
import os
import sys

def plot_spectrometer_with_gif(A_sets, QE_sets, initial_index=0, save_gif=True, gif_name="spectrometer.gif"):
    """
    三軸分光器のXY平面図を描く + GIF 保存機能
    :param A_sets: A1, A2, A3 の角度セットのリスト
    :param QE_sets: スキャン条件のリスト（[ω, h, k, l]）
    :param initial_index: 初期表示する角度セットのインデックス
    :param save_gif: True の場合、GIF を保存する
    :param gif_name: 保存する GIF のファイル名
    """
    # INIファイルの設定読み込み
    config = configparser.ConfigParser()
    if getattr(sys, 'frozen', False):
        ini_path = os.path.join(os.path.dirname(sys.argv[0]), 'config.ini')
    else:
        ini_path = os.path.join(os.path.dirname(__file__), 'config.ini')
    config.read(ini_path)
    
    # 各設定値を取得
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
    max_L = monochromator_to_sample + sample_to_analyzer + analyzer_to_detector

    view_mode = config['settings']['system']

    # プロット準備
    fig, ax = plt.subplots(figsize=(8, 5))
    plt.subplots_adjust(left=0.01, bottom=0.25)

    # スライダー設定
    ax_slider = plt.axes([0.25, 0.10, 0.65, 0.03], facecolor='lightgoldenrodyellow')
    slider = Slider(ax_slider, 'scan number', 1, len(A_sets), valinit=initial_index+1, valstep=1)

    # フレーム保存用リスト
    frames = []

    def update(val):
        index = int(val) - 1
        A1, A2, A3 = A_sets[index]

        x_sample = monochromator_to_sample * np.cos(np.radians(A1))
        y_sample = monochromator_to_sample * np.sin(np.radians(A1))
        x_analyzer = x_sample + sample_to_analyzer * np.cos(np.radians(A1 + A2))
        y_analyzer = y_sample + sample_to_analyzer * np.sin(np.radians(A1 + A2))
        x_detector = x_analyzer + analyzer_to_detector * np.cos(np.radians(A1 + A2 + A3))
        y_detector = y_analyzer + analyzer_to_detector * np.sin(np.radians(A1 + A2 + A3))

        ax.clear()
        ax.scatter(-max_L, 0, color='black', s=10, label='Source')
        ax.add_patch(patches.Circle((0, 0), monochromator_radius, linewidth=2, edgecolor='none', facecolor='blue', linestyle='--'))
        ax.scatter(0, 0, color='blue', s=10, label='Monochromator')
        ax.add_patch(patches.Circle((x_sample, y_sample), sample_goniometer_radius, linewidth=2, edgecolor='none', facecolor='green', linestyle='--'))
        ax.scatter(x_sample, y_sample, color='green', s=10, label='Sample Gonio')
        ax.add_patch(patches.Circle((x_analyzer, y_analyzer), analyzer_radius, linewidth=2, edgecolor='none', facecolor='red', linestyle='--'))
        ax.scatter(x_analyzer, y_analyzer, color='red', s=10, label='Analyzer')
        ax.add_patch(patches.Circle((x_detector, y_detector), detector_radius, linewidth=2, edgecolor='none', facecolor='purple', linestyle='--'))
        ax.scatter(x_detector, y_detector, color='purple', s=10, label='Detector')
        ax.plot([-max_L, 0], [0, 0], color='black', linestyle='--')
        ax.plot([0, x_sample], [0, y_sample], color='blue', linestyle='--')
        ax.plot([x_sample, x_analyzer], [y_sample, y_analyzer], color='green', linestyle='--')
        ax.plot([x_analyzer, x_detector], [y_analyzer, y_detector], color='red', linestyle='--')
        ax.add_patch(Rectangle((floor_position_x, floor_position_y), floor_length, floor_width, linewidth=1, edgecolor='black', facecolor='none', linestyle='--'))

        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim(-max_L, max_L)
        if view_mode == 'left':
            ax.set_ylim(max_L * 1.25, -monochromator_radius)
        else:
            ax.set_ylim(-monochromator_radius, max_L * 1.25)

        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left')

        ax.text(0.4, 1.05, f'ℏω: {QE_sets[index][0]} meV, h: {QE_sets[index][1]}, k: {QE_sets[index][2]}, l: {QE_sets[index][3]}',
                horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

        plt.draw()

        # フレーム保存（GIF用）
        if save_gif:
            fig.canvas.draw()
            image = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
            image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
            frames.append(Image.fromarray(image))

    # スライダー変更時の動作を登録
    slider.on_changed(update)

    # 初期フレームを描画
    update(slider.val)

    if save_gif:
        for val in range(1, len(A_sets) + 1):
            slider.set_val(val)
        frames[0].save(gif_name, save_all=True, append_images=frames[1:], duration=100, loop=0)
        print(f"GIF 保存完了: {gif_name}")

    plt.show()

"""
# サンプルデータ
A_sets = [(0, 30, 60), (10, 40, 70), (20, 50, 80)]
QE_sets = [(10, 1, 0, 0), (15, 1, 1, 0), (20, 1, 1, 1)]

plot_spectrometer_with_gif(A_sets, QE_sets, save_gif=True)
"""