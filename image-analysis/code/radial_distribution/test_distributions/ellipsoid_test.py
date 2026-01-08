"""
Toy example for the radial-position metric using an ellipsoid.

Generates random points inside an ellipsoid and points on its surface, computes
for each point its distance to the centroid and to the closest surface point, and
plots the resulting relative position distribution P = d_centroid / (d_centroid + d_edge).
Also includes a 2D ellipse schematic illustrating the metric.
"""

from qlivecell import compute_distance_xyz_jit, check_or_create_dir, compute_distance_xy_jit, compute_dists_jit
import numpy as np
from numba import njit, prange
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.rcParams.update({
    "text.usetex": True,
})
mpl.rcParams['text.latex.preamble'] = r'\usepackage{siunitx} \sisetup{detect-all} \usepackage{helvet} \usepackage{sansmath} \sansmath'
mpl.rc('font', size=14) 
mpl.rc('axes', labelsize=14) 
mpl.rc('xtick', labelsize=14) 
mpl.rc('ytick', labelsize=14) 
mpl.rc('legend', fontsize=14) 

path_figures = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/radial_distribution/figures/test_examples/"
check_or_create_dir(path_figures)

@njit
def compute_distances(points, points_contour, centroid):
    
    len_points = len(points)
    len_points_contour = len(points_contour)

    dists_contours = np.zeros(len_points)
    dists_centroid = np.zeros(len_points)
    
    for p in prange(len_points):
        point = points[p]
        _dists = np.zeros(len_points_contour)
        for pc in prange(len_points_contour):
            point_c = points_contour[pc]
            dist = compute_distance_xyz_jit(point, point_c)
            _dists[pc] = dist

        dists_contours[p] = np.min(_dists)
        dist = compute_distance_xyz_jit(point, centroid)
        dists_centroid[p] = dist

    return dists_contours, dists_centroid

def generate_points_on_ellipsoid(semi_x=5, semi_y=3, semi_z=2, num_points=100):
    import numpy as np

    # Parameters for the ellipsoid
    a = semi_x  # Semi-major axis length along x-axis
    b = semi_y  # Semi-major axis length along y-axis
    c = semi_z  # Semi-major axis length along z-axis

    # Create a mesh grid for the ellipsoid
    u = np.linspace(0, 2 * np.pi, num_points)
    v = np.linspace(0, np.pi, num_points)
    u, v = np.meshgrid(u, v)

    # Parametric equations for the ellipsoid
    x = a * np.cos(u) * np.sin(v)
    y = b * np.sin(u) * np.sin(v)
    z = c * np.cos(v)

    # Flatten the arrays to create a list of points
    x_flat = x.flatten()
    y_flat = y.flatten()
    z_flat = z.flatten()

    # Combine the points into a single array
    points = np.vstack((x_flat, y_flat, z_flat)).T

    return points


def generate_points_inside_ellipsoid(semi_x=5, semi_y=3, semi_z=2, num_points=10000):
    
    # Parameters for the ellipsoid
    a = semi_x  # Semi-major axis length along x-axis
    b = semi_y  # Semi-major axis length along y-axis
    c = semi_z  # Semi-major axis length along z-axis

    # Generate random points inside a unit sphere using spherical coordinates
    phi = np.random.uniform(0, 2 * np.pi, num_points)
    cos_theta = np.random.uniform(-1, 1, num_points)
    u = np.random.uniform(0, 1, num_points)

    theta = np.arccos(cos_theta)
    r = u ** (1/3)

    # Convert spherical coordinates to Cartesian coordinates
    x_sphere = r * np.sin(theta) * np.cos(phi)
    y_sphere = r * np.sin(theta) * np.sin(phi)
    z_sphere = r * np.cos(theta)

    # Scale points to fit inside the ellipsoid
    x_ellipsoid = a * x_sphere
    y_ellipsoid = b * y_sphere
    z_ellipsoid = c * z_sphere

    # Combine the points into a single array
    points = np.vstack((x_ellipsoid, y_ellipsoid, z_ellipsoid)).T

    return points

def generate_points_on_ellipse(semi_x=5, semi_y=3, num_points=100):
    # Parameters for the ellipse
    a = semi_x  # Semi-major axis length along x-axis
    b = semi_y  # Semi-minor axis length along y-axis

    # Create an array of angles from 0 to 2*pi
    theta = np.linspace(0, 2 * np.pi, num_points)

    # Parametric equations for the ellipse
    x = a * np.cos(theta)
    y = b * np.sin(theta)

    # Combine the points into a single array
    points = np.vstack((x, y)).T

    return points

num_points_base = 1000
radiuses = [[4,1,1]]
DISTS = []
for r, radius in enumerate(radiuses):

    print(radiuses)
    print(radius)
    semi_x, semi_y, semi_z = radius
    num_points = np.rint(num_points_base * np.mean(radius)).astype("int64")

    print(num_points)
    # Example usage:
    points_contour = generate_points_on_ellipsoid(semi_x, semi_y, semi_z, num_points)
    points = generate_points_inside_ellipsoid(semi_x, semi_y, semi_z, num_points)
    xs = points[:,0]
    ys = points[:,1]
    zs = points[:,2]

    centroid = np.mean(points, axis=0)

    dists_contours, dists_centroid = compute_distances(points, points_contour, centroid)

    dists = np.array(dists_centroid) / (np.array(dists_centroid) + np.array(dists_contours))
    DISTS.append(dists)

fig, ax = plt.subplots(1,2, figsize=(8,4))

# Example usage
semi_x = 1
semi_y = 1.5

point_inside = np.array([-0.3, -1.0])
points = generate_points_on_ellipse(semi_x, semi_y)
dists = compute_dists_jit(np.array([point_inside]), np.array(points), compute_distance_xy_jit)
closest = np.argmin(dists, axis=1)
point_contour = points[closest][0]

# Plot the ellipse
ax[0].plot(points[:, 0], points[:, 1], ls='-', label='edge', c='k', lw=4)
ax[0].plot([point_inside[0], point_contour[0]], [point_inside[1], point_contour[1]], c="grey", lw=3)
ax[0].plot([point_inside[0], 0], [point_inside[1], 0], c="grey", lw=3)
ax[0].scatter([point_inside[0]], [point_inside[1]], s=75, label="point", color="brown", zorder=10)
ax[0].scatter([0], [0], s=75, label="centroid", color="k", zorder=10)
ax[0].scatter([point_contour[0]], [point_contour[1]], s=75, label="closest edge", color="yellow", zorder=10)

# Labels and legend
ax[0].title('relative position on ellipse')
ax[0].legend(loc='upper left', bbox_to_anchor=(-0.8, 0.5), frameon=False)
ax[0].set_aspect('equal')
ax[0].spines[['bottom','left', 'right', 'top']].set_visible(False)
ax[0].set_xticks([])
ax[0].set_yticks([])

for r, dists in enumerate(DISTS):
    radius = radiuses[r]
    n_hist, bins, patches = ax[1].hist(dists, alpha=0.5, bins=30, label=r"$R_s$ = {}".format(radius), density=True)
ax[1].legend()
ax[1].set_xlabel(r"relative position on ellipsoid, $P$")
ax[1].spines[['left', 'right', 'top']].set_visible(False)
ax[1].set_yticks([])
plt.tight_layout()
plt.savefig(path_figures+"ellipsoid.svg")
plt.show()


rem=1
bin_n = 50
figg, axx = plt.subplots(1,2)
r = 0
dists = DISTS[r]
radius = radiuses[r]
_counts, bins = np.histogram(dists, bins=bin_n)
axx[0].plot(bins[1:], _counts, ls="-")

dbins = np.mean(np.diff(bins))
bins[1:] -= dbins
bins[0] = 0
counts = [_counts[i]/((4/3)*np.pi*(bins[i+1]**3-bins[i]**3)) for i in range(len(_counts))]
if rem==0:
    axx[1].scatter(bins[1+rem:], counts[rem:], s=30)
else:
    axx[1].scatter(bins[1+rem:-rem], counts[rem:-rem], s=30)
totals = np.sum(_counts)
total_density = np.sum(totals)/((4/3)*np.pi*(bins[-1]**3))
axx[1].plot(bins, np.ones_like(bins)*total_density, color="grey")

plt.show()
plt.scatter(bins[1+rem:-rem], counts[rem:-rem], s=30)
plt.plot(bins, np.ones_like(bins)*total_density, color="grey")
plt.show()