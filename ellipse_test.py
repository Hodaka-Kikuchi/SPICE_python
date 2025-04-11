import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

# 与えられた共分散行列（Qx, Qy, E, Qz の順）
cov = np.array([
    [ 8995.08491475,   409.21439459,   89.09401595,   -0.        ],
    [  409.21439459,  9324.52074736, -1892.89106318,   0.        ],
    [   89.09401595, -1892.89106318,  453.90071286,   -0.        ],
    [   0.,             0.,             0.,          1140.20848813]
])

cov_QxE = cov[[0, 2]][:, [0, 2]]
cov_QyE = cov[[1, 2]][:, [1, 2]]

def plot_cov_ellipse(cov_sub, ax, n_std=1.0, facecolor='none', edgecolor='blue', label=''):
    vals, vecs = np.linalg.eigh(cov_sub)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    width, height = 2 * n_std * np.sqrt(vals)
    ellipse = Ellipse((0, 0), width=width, height=height, angle=theta,
                      facecolor=facecolor, edgecolor=edgecolor, lw=2, label=label)
    ax.add_patch(ellipse)
    return width, height  # サイズを返して表示範囲を調整

# 描画
fig, axs = plt.subplots(1, 2, figsize=(12, 5))

# Qx vs E
axs[0].set_title('Resolution Ellipse in Qx vs E')
w, h = plot_cov_ellipse(cov_QxE, axs[0], edgecolor='blue', label='1σ')
axs[0].axhline(0, color='gray', lw=0.5)
axs[0].axvline(0, color='gray', lw=0.5)
axs[0].set_xlabel('Qx')
axs[0].set_ylabel('E')
axs[0].legend()
axs[0].set_xlim(-w, w)
axs[0].set_ylim(-h, h)
axs[0].set_aspect('auto')  # 'equal' を 'auto' に変更

# Qy vs E
axs[1].set_title('Resolution Ellipse in Qy vs E')
w, h = plot_cov_ellipse(cov_QyE, axs[1], edgecolor='green', label='1σ')
axs[1].axhline(0, color='gray', lw=0.5)
axs[1].axvline(0, color='gray', lw=0.5)
axs[1].set_xlabel('Qy')
axs[1].set_ylabel('E')
axs[1].legend()
axs[1].set_xlim(-w, w)
axs[1].set_ylim(-h, h)
axs[1].set_aspect('auto')  # 同様に auto

plt.tight_layout()
plt.show()
