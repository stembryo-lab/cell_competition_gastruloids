"""
Estimate and correct spectral spillover from a lineage marker (emiRFP or mCherry)
into the p53 channel.

Fits a global linear model C ≈ b0 + s·B using in-cell pixels (per z-plane),
where B is the marker channel and C is p53. The spillover slope s is estimated
globally, while a per-z baseline b0(z) is computed as median(C - s·B).

Saves calibration parameters (s, b0z) for reuse. Corrected p53 per pixel:
C_corr = C - b0(z) - s·B.
"""

BLEED_FROM = "F3"   # <-- set to "F3" or "A12": the channel that bleeds into C ("p53")
C_CHANNEL  = "p53"  # your readout channel (C)

# We set a limit pixels per z-plane when fitting (to keep memory/time bounded).
# Set to None to use ALL in-cell pixels.
SAMPLE_PER_Z = 200_000

from qlivecell import get_file_name, cellSegTrack, check_or_create_dir, get_file_names, fill_channels
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


plt.rcParams.update({"text.usetex": True})
mpl.rcParams['text.latex.preamble'] = r'\usepackage{siunitx} \sisetup{detect-all} \usepackage{helvet} \usepackage{sansmath} \sansmath'
mpl.rc('font', size=18)
mpl.rc('axes', labelsize=18)
mpl.rc('xtick', labelsize=18)
mpl.rc('ytick', labelsize=18)
mpl.rc('legend', fontsize=18)

# ---- Channels (your order) ----
channel_names = ["A12", "p53", "F3", "DAPI"]
ch_B = channel_names.index(BLEED_FROM)
ch_C = channel_names.index(C_CHANNEL)

# =========================
# Helpers
# =========================

def build_union_masks(CT_list):
    """
    Build per-z 2D boolean masks marking in-cell pixels from the union of all cells
    across provided CT objects (e.g., CT_F3 and CT_A12).
    Returns: list of length Z with arrays (Y, X) dtype=bool.
    """
    CT0 = CT_list[0]
    Z = CT0.hyperstack.shape[1]
    Y = CT0.hyperstack.shape[-2]
    X = CT0.hyperstack.shape[-1]
    Mz_list = [np.zeros((Y, X), dtype=bool) for _ in range(Z)]
    for CT in CT_list:
        for cell in CT.jitcells:
            z = int(cell.centers[0][0])
            if z < 0 or z >= Z:
                continue
            # find mask for this z
            try:
                zid = cell.zs[0].index(z)
            except ValueError:
                continue
            mask = cell.masks[0][zid]
            yy = mask[:, 1].astype(np.intp)
            xx = mask[:, 0].astype(np.intp)
            Mz_list[z][yy, xx] = True
    return Mz_list

def sample_indices(n, k):
    if (k is None) or (n <= k):
        return slice(None)
    return np.random.choice(n, size=k, replace=False)

def update_normal_eq_sums(x, y, sums):
    """
    Update sufficient statistics for OLS fit y = b0 + s*x
    sums is a dict with keys: N, sumx, sumy, sumxx, sumxy
    """
    if x.size == 0:
        return
    sums["N"]     += x.size
    sums["sumx"]  += float(np.sum(x))
    sums["sumy"]  += float(np.sum(y))
    sums["sumxx"] += float(np.sum(x * x))
    sums["sumxy"] += float(np.sum(x * y))

def solve_b0_s_from_sums(sums):
    N, sumx, sumy, sumxx, sumxy = sums["N"], sums["sumx"], sums["sumy"], sums["sumxx"], sums["sumxy"]
    denom = (N * sumxx - sumx**2)
    if denom == 0 or N == 0:
        return 0.0, 0.0
    s  = (N * sumxy - sumx * sumy) / denom
    b0 = (sumy - s * sumx) / N
    return float(b0), float(s)

def estimate_b0z_for_file(CT_ref, Mz_list, ch_B, ch_C, s_global):
    """Compute per-z baseline b0(z) = median( C - s_global * B ) over in-cell pixels."""
    Z = CT_ref.hyperstack.shape[1]
    b0z = np.full(Z, np.nan, dtype=np.float64)
    for z in range(Z):
        Mz = Mz_list[z]
        if not np.any(Mz):
            continue
        Bz = CT_ref.hyperstack[0, z, ch_B, :, :].astype(np.float64)
        Cz = CT_ref.hyperstack[0, z, ch_C, :, :].astype(np.float64)
        x  = Bz[Mz].ravel()
        y  = Cz[Mz].ravel()
        resid = y - s_global * x
        b0z[z] = np.median(resid)
    # Fill any empty planes with the median of available planes (fallback to 0)
    if np.any(np.isnan(b0z)):
        if np.any(~np.isnan(b0z)):
            fill = np.nanmedian(b0z)
        else:
            fill = 0.0
        b0z[np.isnan(b0z)] = fill
    return b0z

def correct_cell_pixels(CT_ref, mask, z, ch_B, ch_C, s, b0z):
    """Return per-pixel corrected C for one cell at plane z."""
    yy = mask[:, 1].astype(np.intp)
    xx = mask[:, 0].astype(np.intp)
    C_vals = CT_ref.hyperstack[0, z, ch_C, :, :][yy, xx].astype(np.float32)
    B_vals = CT_ref.hyperstack[0, z, ch_B, :, :][yy, xx].astype(np.float32)
    return C_vals - float(b0z[z]) - float(s) * B_vals

master_path_save = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/p53/segmentation_results/"
path_to_save_figures = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/p53/figures/"
path_to_save_results = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/p53/spillover/"
check_or_create_dir(path_to_save_results)

COND = "KO" 
REP = "n3_Ab"
path_data_dir="/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/p53_analysis/input/{}/{}/".format(COND,REP)
path_save_dir = "{}{}/{}/".format(master_path_save, COND, REP)

files = get_file_names(path_data_dir)
file_to_compute = files[1]

# Accumulate OLS sums across all files/z-planes
sums = {"N": 0, 
        "sumx": 0.0, 
        "sumy": 0.0, 
        "sumxx": 0.0, 
        "sumxy": 0.0
        }

path_data = path_data_dir+file_to_compute
file, embcode = get_file_name(path_data_dir, file_to_compute, allow_file_fragment=False, return_files=False, return_name=True)
path_save = path_save_dir+embcode
check_or_create_dir(path_save)

ch = channel_names.index(BLEED_FROM)
batch_args = {
    'name_format':"ch"+str(ch)+"_{}",
    'extension':".tif",
} 
chans = fill_channels(channel=ch, channel_names=channel_names)

CT = cellSegTrack(
    path_data,
    path_save,
    batch_args=batch_args,
    channels=chans
)
CT.load()
            
# --- Union mask across both populations, per z ---
Mz_list = build_union_masks([CT])

Z = CT.hyperstack.shape[1]
for z in range(Z):
    Mz = Mz_list[z]
    if not np.any(Mz):
        continue
    Bz = CT.hyperstack[0, z, ch_B, :, :].astype(np.float64)
    Cz = CT.hyperstack[0, z, ch_C, :, :].astype(np.float64)

    x = Bz[Mz].ravel()
    y = Cz[Mz].ravel()

    # Optional subsample for speed
    sel = sample_indices(x.size, SAMPLE_PER_Z)
    x = x[sel]; y = y[sel]

    update_normal_eq_sums(x, y, sums)

# Solve for global (session) s and an overall b0 (we'll replace b0 by per-z b0z later)
b0_global, s_global = solve_b0_s_from_sums(sums)
print(f"[Calibration] Estimated global spillover s ({BLEED_FROM} → {C_CHANNEL}): {s_global:.6g}")
print(f"[Calibration] Global intercept b0 (unused; we’ll use per-z): {b0_global:.6g}")

b0z = estimate_b0z_for_file(CT, Mz_list, ch_B, ch_C, s_global)

# Save calibration alongside the file’s outputs (so you can reuse later)
calib_path = os.path.join(path_to_save_results, f"calibration_{BLEED_FROM}_to_{C_CHANNEL}.npz")
np.savez(calib_path, s=s_global, b0z=b0z, bleed_from=BLEED_FROM, c_channel=C_CHANNEL)
print(f"[Saved] {calib_path}  |  s={s_global:.6g}  median(b0z)={np.median(b0z):.6g}")

# Make a small scatter of (C - b0z[z]) vs B for a few planes, before vs after correction
fig, ax = plt.subplots(figsize=(6, 5))
z = 2
Mz = Mz_list[z]
Bz = CT.hyperstack[0, z, ch_B, :, :].astype(np.float64)
Cz = CT.hyperstack[0, z, ch_C, :, :].astype(np.float64)
x = Bz[Mz].ravel()
y = Cz[Mz].ravel()
sel = sample_indices(x.size, 10_000)
x = x[sel]; y = y[sel]
y0 = y - b0z[z]
ycorr = y - b0z[z] - s_global * x
ax.scatter(x, y0, s=10, alpha=0.2, label="before")
ax.scatter(x, ycorr, s=10, alpha=0.2, label="after")
ax.set_xlabel("H2B-emiRFP intensity (B)")
ax.set_ylabel(f"{C_CHANNEL} (C) minus b0(z)")
ax.legend()
plt.tight_layout()
plt.savefig("{}spilloverF3.pdf".format(path_to_save_figures))
plt.show()
