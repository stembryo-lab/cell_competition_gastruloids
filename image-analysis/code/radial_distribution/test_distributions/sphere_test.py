from qlivecell import compute_distance_xyz_jit, check_or_create_dir, compute_distance_xy_jit, compute_dists_jit
import numpy as np
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

def generate_points_on_sphere(num_points, radius=1):
    # Step 1: Generate random points in a cube
    points = np.random.normal(size=(num_points, 3))
    
    # Step 2: Normalize the points to lie on the surface of a unit sphere
    points /= np.linalg.norm(points, axis=1)[:, np.newaxis]
    
    # Step 3: Scale the points to the desired radius
    points *= radius
    
    return points

def generate_points_inside_sphere(num_points, radius=1):
    points = []
    while len(points) < num_points:
        # Generate a random point in a cube
        point = np.random.uniform(-radius, radius, size=3)
        
        # Check if the point is inside the sphere
        if np.linalg.norm(point) <= radius:
            points.append(point)
    
    return np.array(points)

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

num_points_base = 5000
radiuses = [1]
DISTS = []
for r, radius in enumerate(radiuses):
    num_points = num_points_base * radius
    # Example usage:
    points_contour = generate_points_on_sphere(num_points, radius=radius)
    points = generate_points_inside_sphere(num_points, radius=radius)
    xs = points[:,0]
    ys = points[:,1]
    zs = points[:,2]

    centroid = np.mean(points, axis=0)

    dists_contours = []
    dists_centroid = []
    for point in points:
        _dists = []
        for point_c in points_contour:
            dist = compute_distance_xyz_jit(point, point_c)
            _dists.append(dist)

        dists_contours.append(np.min(_dists))
        dist = compute_distance_xyz_jit(point, centroid)
        dists_centroid.append(dist)
        
    dists = np.array(dists_centroid) / (np.array(dists_centroid) + np.array(dists_contours))
    DISTS.append(dists)

# Example usage
semi_x = 1
semi_y = 1

point_inside = np.array([-0.3, -0.6])
points = generate_points_on_ellipse(semi_x, semi_y)
dists = compute_dists_jit(np.array([point_inside]), np.array(points), compute_distance_xy_jit)
closest = np.argmin(dists, axis=1)
point_contour = points[closest][0]

fig, ax = plt.subplots(1,2, figsize=(8,4))
# Plot the ellipse
ax[0].plot(points[:, 0], points[:, 1], ls='-', label='edge', c='k', lw=4)
ax[0].plot([point_inside[0], point_contour[0]], [point_inside[1], point_contour[1]], c="grey", lw=3)
ax[0].plot([point_inside[0], 0], [point_inside[1], 0], c="grey", lw=3)
ax[0].scatter([point_inside[0]], [point_inside[1]], s=75, label="point", color="brown", zorder=10)
ax[0].scatter([0], [0], s=75, label="centroid", color="k", zorder=10)
ax[0].scatter([point_contour[0]], [point_contour[1]], s=75, label="closest edge", color="yellow", zorder=10)

ax[0].title('relative position on ellipse')
ax[0].legend(loc='upper left', bbox_to_anchor=(-0.7, 0.5), frameon=False)
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
plt.savefig(path_figures+"spheres.svg")
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
