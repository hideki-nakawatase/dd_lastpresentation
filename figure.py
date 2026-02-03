import numpy as np
import matplotlib.pyplot as plt
import os


def figure():
    rdir = "el"
    expid = f"{rdir}/ex4-el"

    im, jm, km = 100, 62, 12
    dx, dy, dz = 200.0e3, 100.0e3, 50.0

    nbgn, nend, nskp = 36, 36, 1

    if not os.path.exists(rdir):
        os.makedirs(rdir)

    # 構造体定義（Fortranバイナリ読込用）
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

        time = chunk["time"][0]
        # order='F' で読み込み
        u = chunk["u"][0].reshape((im + 2, jm + 2, km + 2), order="F")
        v = chunk["v"][0].reshape((im + 2, jm + 2, km + 2), order="F")
        e = chunk["e"][0].reshape((im + 2, jm + 2), order="F")

        # 座標の定義
        xc = (np.arange(im + 2) - 0.5) * dx / 1000.0
        yc = (np.arange(jm + 2) - 0.5) * dy / 1000.0
        zc = (np.arange(km + 2) - 0.5 - km) * dz

        # --- 1. 流速ベクトル図 (XY断面) ---
        target_k = km - 2
        # C-格子からセルセンターへ（内点のみ）
        uc = 0.5 * (u[0:-1, 1:-1, target_k] + u[1:, 1:-1, target_k]) * 100
        vc = 0.5 * (v[1:-1, 0:-1, target_k] + v[1:-1, 1:, target_k]) * 100
        # サイズを合わせるために調整
        uc = uc[0:im, 0:jm]
        vc = vc[0:im, 0:jm]
        u_abs = np.sqrt(uc**2 + vc**2)

        XC, YC = np.meshgrid(xc[1:-1], yc[1:-1], indexing="ij")

        fig1, ax1 = plt.subplots(figsize=(8, 6))
        ax1.contour(xc, yc, e.T, colors="gray", alpha=0.5)
        # データの中心付近を表示（格子サイズに合わせてスライス）
        Q = ax1.quiver(XC, YC, uc, vc, u_abs, cmap="jet", pivot="mid")
        plt.colorbar(Q, label="Velocity [cm/s]", shrink=0.8)
        ax1.set_title(f"Velocity Field (Time: {time:.1f}s)")
        ax1.set_xlabel("X [km]")
        ax1.set_ylabel("Y [km]")
        fig1.savefig(f"{rdir}/velocity_xy_{n:06}.png", dpi=200)
        plt.close(fig1)

        # --- 2. 流関数 (Stream function) ---
        # 鉛直積分流量 (内点のみ)
        um = np.sum(u[1 : im + 1, 1 : jm + 1, 1 : km + 1], axis=2) * dz
        psi = np.zeros((im, jm))
        # y方向へ積分
        psi[:, 1:] = np.cumsum(um[:, :-1] * dy, axis=1)

        xc_psi = xc[1:-1]
        yc_psi = yc[1:-1]
        X_P, Y_P = np.meshgrid(xc_psi, yc_psi, indexing="ij")

        fig2, ax2 = plt.subplots(figsize=(8, 6))
        cntr = ax2.contour(X_P, Y_P, psi, levels=15, cmap="RdBu_r")
        ax2.clabel(cntr, fmt="%.1e")
        ax2.set_title("Stream Function")
        ax2.set_aspect("equal")
        fig2.savefig(f"{rdir}/stream_func_{n:06}.png", dpi=200)
        plt.close(fig2)

        # --- 3. 鉛直断面図 (YZ断面 at x=im/2) ---
        # セルセンターでのu成分
        uc_yz = 0.5 * (u[im // 2, :, :] + u[im // 2 + 1, :, :]) * 100
        # 座標メッシュ作成
        YC_YZ, ZC_YZ = np.meshgrid(yc, zc, indexing="ij")

        fig3, ax3 = plt.subplots(figsize=(8, 6))
        # uc_yz は (jm+2, km+2) のサイズ
        pc = ax3.pcolormesh(
            YC_YZ, ZC_YZ, uc_yz, cmap="RdBu_r", shading="auto", vmin=-70, vmax=70
        )
        plt.colorbar(pc, label="u [cm/s]")
        ax3.set_title(f"YZ Section of U-velocity (x={xc[im//2]:.0f}km)")
        ax3.set_xlabel("y [km]")
        ax3.set_ylabel("z [m]")
        ax3.set_ylim(-600, 0)
        fig3.savefig(f"{rdir}/figure_yz_{n:06}.png")
        plt.close(fig3)

    # --- 4 風応力の分布 ---
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

    fig4, ax4 = plt.subplots(figsize=(8, 6))
    Qw = ax4.quiver(XW, YW, tx, ty, pivot="mid", color="blue")
    ax4.set_title("Wind Stress Distribution [N/m^2]")
    ax4.set_xlabel("x [km]")
    ax4.set_ylabel("y [km]")
    fig4.savefig(f"{rdir}/wind_stress.png")
    plt.close(fig4)


if __name__ == "__main__":
    figure()
