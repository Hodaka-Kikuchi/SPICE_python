import numpy as np

def resolution_convolution(M, H, K, L, W,
                           tqx, tqyy, tqyx, tqwx, tqwy, tqww, tqz,
                           xvec, yvec, zvec,
                           sqw, p,
                           prefactor, detM):
    """
    MATLABコードを移植した分解能畳み込み（Python/NumPy版）。

    入力（想定）：
      M         : tuple/list (M1, M2)  (整数)
      H,K,L,W   : 1D arrays, 長さ N (各スキャン点)
      tqx,tqyy,tqyx,tqwx,tqwy,tqww,tqz : 1D arrays, 長さ N (tq係数)
      xvec,yvec,zvec : arrays shape (3, N) または (N,3)（下で対応）
                     MATLABコードは xvec(1,i) 等でアクセスする前提
                     この実装では xvec[0,i], xvec[1,i], xvec[2,i] を使用
      sqw       : callable: sqw(H1, K1, L1, W1, p) -> array shape (modes, npts)
                  （npts はグリッド数 = (2*M1+1)^3）
      p         : パラメータ（sqw に渡す）
      prefactor : array shape (modes, N)
      detM      : 行列 M の行列式（スカラー）

    出力：
      conv : 1D array 長さ N（各スキャン点の畳み込み結果）
    """

    M1, M2 = M[0], M[1]
    # steps and grids (MATLAB と同じ定義)
    step1 = np.pi / (2 * M1 + 1)
    step2 = np.pi / (2 * M2 + 1)
    dd1 = np.linspace(-np.pi/2 + step1/2, np.pi/2 - step1/2, 2*M1 + 1)
    dd2 = np.linspace(-np.pi/2 + step2/2, np.pi/2 - step2/2, 2*M2 + 1)

    modes = prefactor.shape[0]
    N = len(H)

    # ndgrid(dd1,dd1,dd1) 相当 (i,j,k) -> shape (n1,n1,n1)
    n1 = dd1.size
    cx, cy, cw = np.meshgrid(dd1, dd1, dd1, indexing='ij')  # shape (n1,n1,n1)

    # flatten to 1D arrays of length n1**3 (MATLAB の cx(1:end) 相当)
    tx = np.tan(cx.ravel())
    ty = np.tan(cy.ravel())
    tw = np.tan(cw.ravel())

    tz = np.tan(dd2)  # shape (n2,)

    # 正規化係数（MATLAB と同じ式）
    norm = np.exp(-0.5 * (tx**2 + ty**2)) * (1 + tx**2) * (1 + ty**2) \
           * np.exp(-0.5 * (tw**2)) * (1 + tw**2)   # shape (n1**3,)
    normz = np.exp(-0.5 * (tz**2)) * (1 + tz**2)     # shape (n2,)

    # 結果格納
    convs = np.zeros((modes, N), dtype=float)
    conv = np.zeros(N, dtype=float)

    # ループは MATLAB と同じ構造 (iz -> dd2 のループ, i -> スキャン点ループ)
    # xvec,yvec,zvec は shape (3,N) を想定（もし (N,3) の場合は転置して渡すこと）
    x1 = xvec[0, :]  # shape (N,)
    x2 = xvec[1, :]
    x3 = xvec[2, :]

    y1 = yvec[0, :]
    y2 = yvec[1, :]
    y3 = yvec[2, :]

    z1 = zvec[0, :]
    z2 = zvec[1, :]
    z3 = zvec[2, :]

    # 主要ループ
    npoints = tx.size  # = n1**3
    for iz_idx, tz_val in enumerate(tz):
        nz_factor = normz[iz_idx]   # scalar for this iz

        for i in range(N):
            # 各 t* ベクトルは長さ npoints の配列
            dQ1 = tqx[i] * tx                      # shape (npoints,)
            dQ2 = tqyy[i] * ty + tqyx[i] * tx
            dW  = tqwx[i] * tx + tqwy[i] * ty + tqww[i] * tw
            dQ4 = tqz[i] * tz_val                  # scalar

            # H1,K1,L1,W1 は長さ npoints の配列
            H1 = H[i] + dQ1 * x1[i] + dQ2 * y1[i] + dQ4 * z1[i]
            K1 = K[i] + dQ1 * x2[i] + dQ2 * y2[i] + dQ4 * z2[i]
            L1 = L[i] + dQ1 * x3[i] + dQ2 * y3[i] + dQ4 * z3[i]
            W1 = W[i] + dW

            # sqw が (modes, npoints) を返すことを想定
            # 例: out = sqw(H1, K1, L1, W1, p)
            out = sqw(H1, K1, L1, W1, p)   # shape (modes, npoints)

            # 各モードごとに重みを足しこむ
            # add = out[j,:] * norm * nz_factor
            # convs[j,i] += sum(add)
            # これを vectorized に（モード方向ループ）
            # out * norm  -> shape (modes, npoints) broadcasting
            weighted = out * norm  # broadcasting: norm (npoints,) applied to axis 1
            # sum over npoints
            sums = weighted.sum(axis=1) * nz_factor  # shape (modes,)
            convs[:, i] += sums

        # end for i
    # end for iz

    # conv(i) = sum(convs(:,i).*prefactor(:,i))
    conv = np.sum(convs * prefactor, axis=0)   # shape (N,)

    # 最後のスケール
    conv = conv * (step1**3) * step2 / np.sqrt(detM)

    return conv
