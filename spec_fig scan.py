import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
from matplotlib.widgets import Slider

def plot_spectrometer(A1, A2, A3):
    """
    三軸分光器のXY平面図を描く
    :param A1: モノクロメータとサンプルゴニオのなす角度 (度)
    :param A2: サンプルゴニオとアナライザのなす角度 (度)
    :param A3: アナライザとディテクターのなす角度 (度)
    """
    A1 = A1
    A2 = -A2
    A3 = A3
    
    # 距離と半径の設定
    monochromator_radius = 800  # モノクロメータの半径
    monochromator_to_sample = 1600  # モノクロメータとサンプルゴニオ間距離
    sample_goniometer_radius = 350  # サンプルゴニオの半径
    sample_to_analyzer = 700  # サンプルゴニオとアナライザの間距離
    analyzer_radius = 250  # アナライザの半径
    analyzer_to_detector = 500  # アナライザとディテクターの間距離
    detector_radius = 300  # ディテクターの半径
    
    # 距離のマックス
    max_L = monochromator_to_sample + sample_to_analyzer + analyzer_to_detector
    
    # モノクロメータの位置 (原点)
    x_mono = 0
    y_mono = 0

    # プロット設定
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)
    
    # ソース
    ax.scatter(-max_L, 0, color='black', s=1, label='Source')
    
    # スライダー設定
    ax_slider = plt.axes([0.25, 0.01, 0.65, 0.03], facecolor='lightgoldenrodyellow')
    slider = Slider(ax_slider, 'Angle Set', 0, len(A_sets)-1, valinit=0, valstep=1)

    def update(val):
        # スライダーの値に基づいて角度セットを選択
        angle_set = A_sets[int(val)]
        A1, A2, A3 = angle_set
        
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
        
        # floorの外形を書く。
        ax.add_patch(Rectangle((-1850, 800), 4330, 2375, linewidth=1, edgecolor='black', facecolor='none', linestyle='--'))
        
        # 軸設定
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim(-max_L, max_L)
        ax.set_ylim(-800, max_L*1.25)

        ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
        plt.draw()

    slider.on_changed(update)

    # 初期の描画を行う
    update(slider.val)
    
    def on_key(event):
        """
        左右矢印キーでスライダーを動かす
        """
        if event.key == 'right':
            slider.set_val(min(slider.val + 1, len(A_sets) - 1))
        elif event.key == 'left':
            slider.set_val(max(slider.val - 1, 0))
    
    # キーイベントを設定
    fig.canvas.mpl_connect('key_press_event', on_key)

    plt.show()

# 角度のセットを定義
A_sets = [
    [50, -90, 90],  # セット 1
    [55, -80, 90],  # セット 2
    [60, -70, 90],  # セット 3
    [65, -60, 90],  # セット 4
    [70, -50, 90],  # セット 5
    [75, -40, 90],  # セット 6
    [80, -30, 90],  # セット 7
    [85, -20, 90],  # セット 8
    [90, -10, 90],  # セット 9
] 

# 初期の角度設定
plot_spectrometer(50, 30, 90)  # 初期の描画を行う
