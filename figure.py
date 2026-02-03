import numpy as np
import matplotlib.pyplot as plt
import os


def figure():
    rdir = "ty_applied"
    expid = f"{rdir}/ex4-ty"

    # 計算パラメータの統一（km, jmなどがセクションごとに変わる場合は注意が必要ですが、ここでは統合します）
    im, jm, km = 100, 62, 12
    dx, dy, dz = 200.0e3, 100.0e3, 50.0  # 基本のdzを使用

    nbgn, nend, nskp = 36, 36, 1

    # 保存先ディレクトリの作成
    if not os.path.exists(rdir):
        os.makedirs(rdir)

    # 共通のデータ型定義
    dtyp = np.dtype(
        [
            ("time", "<d"),
            ("u", f"<{(im+2)*(jm+2)*(km+2)}d"),
            ("v", f"<{(im+2)*(jm+2)*(km+2)}d"),
            ("e", f"<{(im+2)*(jm+2)}d"),
            ("w", f"<{(im+2)*(jm+2)*(km+2)}d"),
            ("p", f"<{(im+2)*(jm+2)*(km+2)}d"),
            ("t", f"<{(im+2)*(jm+2)*(km+2)}d"),
            ("s", f"<{(im+2)*(jm+2)*(km+2)}d"),
            ("r", f"<{(im+2)*(jm+2)*(km+2)}d"),
        ]
    )

    for n in range(nbgn, nend + 1, nskp):
        fname = f"{expid}.n{n:06}"
        try:
            with open(fname, "rb") as fp:
                chunk = np.fromfile(fp, dtype=dtyp)
        except FileNotFoundError:
            print(f"Skipping: {fname} not found.")
            continue

        # データ展開 (order='F' for Fortran-like binary)
        time = chunk["time"][0]
        u = chunk["u"][0].reshape((im + 2, jm + 2, km + 2), order="F")
        v = chunk["v"][0].reshape((im + 2, jm + 2, km + 2), order="F")
        e = chunk["e"][0].reshape((im + 2, jm + 2), order="F")

        # --- 1. 流速ベクトル図 (XY断面) ---
        # Staggered格子からセルセンターへの変換 (スライスによる高速化)
        uc_raw = np.zeros((im + 2, jm + 2, km + 2))
        vc_raw = np.zeros((im + 2, jm + 2, km + 2))
        uc_raw[:, 1 : jm + 1, :] = 0.5 * (u[:, 0:jm, :] + u[:, 1 : jm + 1, :])
        vc_raw[:, 1 : jm + 1, :] = 0.5 * (v[:, 0:jm, :] + v[:, 1 : jm + 1, :])

        # 特定の深さ (km-2) を抽出して単位変換
        target_k = km - 2
        uc = uc_raw[1:-1, 1:-1, target_k] * 100
        vc = vc_raw[1:-1, 1:-1, target_k] * 100
        u_abs = np.sqrt(uc**2 + vc**2)

        xc = (np.arange(im + 2) - 0.5) * dx / 1000.0
        yc = (np.arange(jm + 2) - 0.5) * dy / 1000.0
        XC, YC = np.meshgrid(xc[1:-1], yc[1:-1], indexing="ij")

        fig1, ax1 = plt.subplots(figsize=(8, 6))
        # e (海面水位など) のコンター
        ax1.contour(xc, yc, e.T, colors="gray", alpha=0.5)  # eのサイズに合わせる
        Q = ax1.quiver(XC, YC, uc, vc, u_abs, cmap="jet", pivot="mid")
        plt.colorbar(Q, label="Velocity [cm/s]", shrink=0.8)
        ax1.set_title(f"Velocity Field (Time: {time:.1f}s)")
        ax1.set_xlabel("X [km]")
        ax1.set_ylabel("Y [km]")
        ax1.set_aspect("auto")
        fig1.savefig(f"{rdir}/velocity_xy_{n:06}.png", dpi=200)
        plt.close(fig1)

        # --- 2. 流関数 (Stream function) ---
        # 鉛直積分流量
        um = np.sum(u[1 : im + 1, 1 : jm + 1, 1 : km + 1], axis=2) * dz
        psi = np.zeros((im, jm))
        psi[:, 1:] = np.cumsum(um[:, :-1] * dy, axis=1)  # ループをcumsumに置換

        xc_psi = np.arange(im) * dx / 1000.0
        yc_psi = np.arange(jm) * dy / 1000.0
        X_P, Y_P = np.meshgrid(xc_psi, yc_psi, indexing="ij")

        fig2, ax2 = plt.subplots(figsize=(8, 6))
        cntr = ax2.contour(X_P, Y_P, psi, levels=15, cmap="RdBu_r")
        ax2.clabel(cntr, fmt="%.1e")
        ax2.set_title("Stream Function")
        ax2.set_aspect("equal")
        fig2.savefig(f"{rdir}/stream_func_{n:06}.png", dpi=200)
        plt.close(fig2)

    # --- 3. 風応力の分布 (定数分布なのでループ外) ---
    im_w, jm_w = 20, 13
    t0 = -0.1
    tx = np.zeros((im_w, jm_w))
    ty = np.zeros((im_w, jm_w))
    xw = np.linspace(0, im * dx / 1000.0, im_w)
    yw = np.linspace(0, jm * dy / 1000.0, jm_w)
    XW, YW = np.meshgrid(xw, yw, indexing="ij")

    # 風応力の計算（サイン・コサイン分布）
    for j in range(jm_w):
        y_rel = (j - jm_w / 2) / (jm_w / 2)
        tx[:, j] = t0 * np.cos(np.pi * y_rel)
        ty[:, j] = t0 * np.sin(np.pi * y_rel / 2)

    fig3, ax3 = plt.subplots(figsize=(8, 6))
    Qw = ax3.quiver(XW, YW, tx, ty, pivot="mid", color="blue")
    ax3.set_title("Wind Stress Distribution [N/m^2]")
    ax3.set_xlabel("x [km]")
    ax3.set_ylabel("y [km]")
    fig3.savefig(f"{rdir}/wind_stress.png")
    plt.close(fig3)


if __name__ == "__main__":
    figure()
