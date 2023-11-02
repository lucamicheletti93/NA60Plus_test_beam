import random
from scipy import optimize

# Generate a list of 10 random integers between 1 and 100
ntrk = 100
mu_vx = 2
mu_vy = 3
mu_px = -1
mu_py = 0

sigma_vx = 2
sigma_vy = 3
sigma_px = 1
sigma_py = 1
#reconstructed tracks
vx_l = [random.gauss(mu_vx, sigma_vx) for _ in range(ntrk)]
vy_l = [random.gauss(mu_vy, sigma_vy) for _ in range(ntrk)]
px_l = [random.gauss(mu_px, sigma_px) for _ in range(ntrk)]
py_l = [random.gauss(mu_py, sigma_py) for _ in range(ntrk)]
#initial guess from the reconstructed beam
xpb = 2
ypb = 2

def distToPvSquared(parameters):
    xv, yv, zv = parameters
    dist = 0
    for vx,vy,px,py in zip(vx_l, vy_l, px_l, py_l):
        dist2 += (vx*zv+px-xv)**2+(vy*zv+py-yv)**2
    return dist

# Initial guess
initial_guess = [xpb, ypb ,0]

# Optimize the function using the BFGS method
result = optimize.minimize(objective_function, initial_guess, method='BFGS')

# Retrieve the results
optimized_parameters = result.x
minimum_value = result.fun
print(optimized_parameters)
