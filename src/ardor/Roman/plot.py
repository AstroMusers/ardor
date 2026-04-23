#%%
from matplotlib import pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u
from astropy.table import Table
from astropy.visualization import astropy_mpl_style
import matplotlib.font_manager as fm
import matplotlib.patheffects as pe
def lm_font():
    font_path = '/ugrad/whitsett.n/fonts/latin-modern-roman/lmroman10-regular.otf'
    fm.fontManager.addfont(font_path)
    plt.rcParams['font.family'] = 'Latin Modern Roman'
    plt.rcParams['font.sans-serif'] = ['Latin Modern Roman']
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Latin Modern Roman'
lm_font()

# # Field geometry (Galactic)
# l0 = 0.5 * u.deg
# b0 = -1.4 * u.deg
# half_l = (23 * 5 * u.arcmin) / 2
# half_b = (45 * u.arcmin) / 2

# lmin = (l0 - half_l).to(u.deg)
# lmax = (l0 + half_l).to(u.deg)

# plt.style.use(astropy_mpl_style)

# # Load CSV
# csv_path = "/ugrad/whitsett.n/Roman/gaia_results.csv"
# tbl = Table.read(csv_path, format="csv")

# # ------------------------------------------------------------------
# # 1) Convert RA/Dec -> l/b and save back to CSV
# # ------------------------------------------------------------------
# if ("l" not in tbl.colnames) or ("b" not in tbl.colnames):
#     ra_all = np.array(tbl["ra"], dtype=float)
#     dec_all = np.array(tbl["dec"], dtype=float)

#     good = np.isfinite(ra_all) & np.isfinite(dec_all)
#     l_all = np.full(len(tbl), np.nan, dtype=float)
#     b_all = np.full(len(tbl), np.nan, dtype=float)

#     gal = SkyCoord(ra=ra_all[good] * u.deg, dec=dec_all[good] * u.deg, frame="icrs").galactic
#     l_all[good] = gal.l.deg
#     b_all[good] = gal.b.deg

#     tbl["l"] = l_all
#     tbl["b"] = b_all
#     tbl.write(csv_path, format="csv", overwrite=True)
#     print(f"Updated CSV with l/b columns: {csv_path}")

# # ------------------------------------------------------------------
# # 2) Plot in Galactic coordinates and filter to S1..S5
# # ------------------------------------------------------------------
# l = np.array(tbl["l"], dtype=float)
# b = np.array(tbl["b"], dtype=float)
# teff = np.array(tbl["teff_gspphot"], dtype=float)

# # Wrap l for continuity around 0/360
# l_wrap = Angle(l * u.deg).wrap_at(180 * u.deg).deg
# lmax_w = Angle(lmax).wrap_at(180 * u.deg).deg
# b0_deg = b0.to_value(u.deg)

# # Base filter
# m = np.isfinite(l_wrap) & np.isfinite(b) & np.isfinite(teff) & (teff < 4000)
# l_wrap, b, teff = l_wrap[m], b[m], teff[m]

# # Build 5 adjacent fields (23' x 45') and union mask
# # Build 5 adjacent fields (23' x 45')
# field_w = (23 * u.arcmin).to_value(u.deg)   # along l
# field_h = (45 * u.arcmin).to_value(u.deg)   # along b
# half_h = field_h / 2.0

# field_bounds = []  # (l_lo, l_hi, b_lo, b_hi)
# for i in range(5):
#     l_hi = lmax_w - i * field_w
#     l_lo = lmax_w - (i + 1) * field_w
#     b_lo = b0_deg - half_h
#     b_hi = b0_deg + half_h
#     field_bounds.append((l_lo, l_hi, b_lo, b_hi))

# # Roman-like detector layout inside each field: 3 x 6 = 18 SCAs
# nx, ny = 3, 6
# gap_arcsec = 0.0   # set e.g. 10.0 to include small inter-detector gaps
# gap_deg = gap_arcsec / 3600.0

# detector_bounds = []  # (field_id, det_id, dl_lo, dl_hi, db_lo, db_hi)
# in_any_detector = np.zeros_like(l_wrap, dtype=bool)

# for fidx, (l_lo, l_hi, b_lo, b_hi) in enumerate(field_bounds, start=1):
#     det_w = (field_w - (nx - 1) * gap_deg) / nx
#     det_h = (field_h - (ny - 1) * gap_deg) / ny

#     det_id = 1
#     for ix in range(nx):
#         dl_lo = l_lo + ix * (det_w + gap_deg)
#         dl_hi = dl_lo + det_w
#         for iy in range(ny):
#             db_lo = b_lo + iy * (det_h + gap_deg)
#             db_hi = db_lo + det_h

#             detector_bounds.append((fidx, det_id, dl_lo, dl_hi, db_lo, db_hi))
#             in_any_detector |= (
#                 (l_wrap >= dl_lo) & (l_wrap <= dl_hi) &
#                 (b >= db_lo) & (b <= db_hi)
#             )
#             det_id += 1

# # Keep only stars in detector footprints
# l_wrap, b, teff = l_wrap[in_any_detector], b[in_any_detector], teff[in_any_detector]


# teff_plot = np.clip(teff, None, 6000.0)

# fig, ax = plt.subplots(figsize=(7, 4))
# sc = ax.scatter(
#     l_wrap, b,
#     c=teff_plot,
#     s=0.1,
#     alpha=0.7,
#     cmap="plasma_r",
#     vmin=float(np.min(teff_plot)),
#     vmax=6000.0
# )

# # Draw detector outlines + field outlines/labels
# for fidx, det_id, dl_lo, dl_hi, db_lo, db_hi in detector_bounds:
#     ax.plot(
#         [dl_lo, dl_hi, dl_hi, dl_lo, dl_lo],
#         [db_lo, db_lo, db_hi, db_hi, db_lo],
#         color="white",
#         lw=0.45,
#         alpha=0.35,
#         zorder=2
#     )

# for i, (l_lo, l_hi, b_lo, b_hi) in enumerate(field_bounds, start=1):
#     ax.plot(
#         [l_lo, l_hi, l_hi, l_lo, l_lo],
#         [b_lo, b_lo, b_hi, b_hi, b_lo],
#         color="deepskyblue",
#         lw=2.0,
#         alpha=0.95,
#         zorder=3
#     )
#     txt = ax.text(
#         (l_lo + l_hi) / 2.0, b0_deg, f"S{i}",
#         fontsize=13,
#         fontweight="bold",
#         color="yellow",
#         ha="center",
#         va="center",
#         zorder=4
#     )
#     txt.set_path_effects([pe.withStroke(linewidth=2.5, foreground="black")])

# ax.set_xlabel("Galactic l [deg]")
# ax.set_ylabel("Galactic b [deg]")
# ax.grid(True, alpha=0.3)

# # Bounds of plot = bounds of plotted data
# if l_wrap.size > 0 and b.size > 0:
#     # keep RA-like visual convention (decreasing to the right)
#     ax.set_xlim(np.max(l_wrap), np.min(l_wrap))
#     ax.set_ylim(np.min(b), np.max(b))

# # Optional: keep visual convention similar to RA plots
# ax.invert_xaxis()

# cbar = plt.colorbar(sc, ax=ax)
# cbar.set_label("Teff [K]")

# plt.tight_layout()
# plt.savefig("/ugrad/whitsett.n/Roman/gaia_stars_s1_s5.png", dpi=400)
# plt.show()