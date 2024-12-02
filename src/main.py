# Zero-pressure gradient (ZPG) boundary layer
#%%
# Calculate the velocity and shear stress distribution
# Import necessary modules
import numpy as np

def ascii_loading_bar(j, N, bar_length=50):
    # Calculate the percentage and number of filled characters
    percent = (j / N) * 100
    filled_length = int(bar_length * j // N)  # How much of the bar should be filled
    bar = '=' * filled_length + '-' * (bar_length - filled_length)  # Build the bar
    
    # Print the loading bar with percentage
    print(f"\r[{bar}] {percent:.2f}%", end='', flush=True)

print(
    f"Initializing the simulation process...\n"
    f"Configuring parameters and preparing to solve the system.\n"
)

# Define parameters
Re = 10**3  # Reynolds number
L = 10      # Length in x direction
H = 2       # Height in y direction
N = 100     # Number of points in x direction
M = 600     # Number of points in y direction

# Initialize the grid:
# N + 1 points in the x direction (including initial value @ x = 0)
# M + 1 points in the y direction (including boundary @ y = 0)

# The interior points span a matrix of (M)x(N)

# Define grid spacing 
x = np.linspace(0, L, N + 1)
dx = x[1]-x[0]
y = np.linspace(0, H, M + 1)
dy = y[1]-y[0]

# Define initial conditions functions
def get_u0(y):
    return np.tanh(5*y)

def get_v0(y):
    return 0

def get_t0(y):
    return 5 / (Re * np.cosh(5 * y)**2)

# Define the solution matrices
U = np.zeros([M + 1, N + 1])
V = np.zeros([M + 1, N + 1])
T = np.zeros([M + 1, N + 1])

# Calculate initial condition
U[:,0] = get_u0(y)
V[:,0] = get_v0(y)
T[:,0] = get_t0(y)

# Define constants for the A matrix
a1 = 1 / (2 * dx)
a2 = -1 / (Re * dy)
b1 = 1 / (2 * dx)
b2 = 1 / (Re * dy)
c1 = 1 / (2 * dy)
d1 = -1 / (2 * dy)
e2 = 1 / 2
e3 = -1 / (2 * dy)
f2 = 1 / 2
f3 = 1 / (2 * dy)

# Build the static part of the A matrix
A = np.zeros([3 * M, 3 * M])

# Loop to populate the interior points of matrix A
for k in range(1, M - 1):
    # A1 block
    A[k, k - 1:k + 1] = [b1, a1]

    # A2 block
    A[k + M, k - 1:k + 1] = [b2, a2]

    # B1 block
    A[k, (k - 1 + M - 1):(k + 1 + M - 1)] = [d1, c1]

    # C2 block
    A[k + M, (k + 2 * M - 1):(k + 2 + 2 * M - 1)] = [f2, e2]

    # C3 block
    A[k + 2 * M, (k + 2 * M - 1):(k + 2 + 2 * M - 1)] = [f3, e3]

# Boundary equations
A[0, 0] = a1
A[M, 0] = a2
A[M-1, M-2] = b1
A[2*M-1, M-2] = b2
A[0, M-1] = c1
A[M-1, (2*M-3):(2*M-1)] = [d1, c1]
A[M, (2*M-1):(2*M+1)] = [f2, e2]
A[2 * M - 1, -2:] = [f2, e2]
A[2*M, (2*M-1):(2*M+1)] = [f3, e3]
A[3 * M - 1, -2:] = [f3, e3]

# Initialize vectors for the iterative solution
Ui = np.zeros(M + 1)
Vi = np.zeros(M + 1)
Ti = np.zeros(M + 1)

# RHS
g = np.zeros(3 * M)
# Solution vector (dU, dV, dT)
X = np.zeros(3 * M)

iterate = 1
it = np.zeros(N + 1)    # vector which saves the iteration count per step

print(
    f"System setup complete. Starting the iterative solver...\nThis process may take some time.\n"
)

# Calculate values along x direction
for j in range(1, N + 1):
    ascii_loading_bar(j, N)  # Update the loading bar with the current progress  

    # Initial iterate - solution from previous line
    Ui[:] = U[:, j - 1]
    Vi[:] = V[:, j - 1]
    Ti[:] = T[:, j - 1]

    iterate = 1
    while iterate > 0:  # Iterate to refine the solution  
        # Build the right-hand side (RHS)
        # g1[:]
        g[:M] = (a1 * U[1:, j - 1] + b1 * U[:-1, j - 1] - c1 * V[1:, j - 1] - d1 * V[:-1, j - 1] 
                - a1 * Ui[1:] - b1 * Ui[:-1] - c1 * Vi[1:] - d1 * Vi[:-1])
        
        # g2[:]
        g[M:2*M] = -e2 * Ti[1:] - f2 * Ti[:-1] - a2 * Ui[1:] - b2 * Ui[:-1]

        # g3[:]
        g[2*M:] = (1 / (2 * dx) * (U[1:, j - 1]**2 + U[:-1, j - 1]**2) 
                - 1 / (2 * dy) * (U[1:, j - 1] * V[1:, j - 1] - U[:-1, j - 1] * V[:-1, j - 1]) 
                + 1 / (2 * dy) * (T[1:, j - 1] - T[:-1, j - 1])
                - 1 / (2 * dx) * (Ui[1:]**2 + Ui[:-1]**2)
                - 1 / (2 * dy) * (Ui[1:] * Vi[1:] - Ui[:-1] * Vi[:-1]) 
                + 1 / (2 * dy) * (Ti[1:] - Ti[:-1]))

        # Update coefficients of the A matrix 
        for k in range(1, M - 1):
            a3 = 1 / dx * Ui[k+1] + 1 / (2 * dy) * Vi[k+1]
            b3 = 1 / dx * Ui[k] - 1 / (2 * dy) * Vi[k]
            c3 = 1 / (2 * dy) * Ui[k+1]
            d3 = -1 / (2 * dy) * Ui[k]

            # A3 block
            A[k + 2 * M, k - 1 : k + 1] = [b3, a3]
            # B3 block
            A[k + 2 * M, (k - 1 + M - 1):(k + 1 + M - 1)] = [d3, c3]
        
        # Boundary equations for A3/B3 equations
        A[2*M, 0] = 1 / dx * Ui[1] + 1 / (2 * dy) * Vi[1]
        A[3*M-1, M-2] = 1 / dx * Ui[M - 1] - 1 / (2 * dy) * Vi[M - 1]
        A[2*M,M-1] = 1 / (2 * dy) * Ui[1]
        A[3*M-1, (2*M-3):(2*M-1)] = [-1 / (2 * dy) * Ui[M - 1], 1 / (2 * dy) * Ui[M - 1]]
        
        # Solve the linear system
        X = np.linalg.solve(A, g)

        # Update the iterate values
        # U(i+1) = U(i) + dU        (U(0)=0, U(M)=Ue)
        Ui[1:-1] = Ui[1:-1] + X[:M-1]
        Ui[0] = 0
        Ui[-1] = 1
        # V(i+1) = V(i) + dV        (V(0)=0)
        Vi[1:] = Vi[1:] + X[M-1:2*M-1]
        Vi[0] = 0
        # T(i+1) = T(i) + dT
        Ti[:] = Ti[:] + X[2*M-1:]
        
        it[j] += 1
        # Check if max iterative value of u, v and t are below some threshold values
        if (np.amax(X[:M-1]) < 10**(-3) and
            np.amax(X[M-1:2*M-1]) < 10**(-3) and
            np.amax(X[2*M-1:]) < 10**(-6)):
            iterate = 0
    # Save the calculated values into the solution matrix
    U[:, j] = Ui[:]
    V[:, j] = Vi[:]
    T[:, j] = Ti[:]

# Boundary layer thickness
d99 = np.zeros(N + 1)
# Displacement thickness
d = np.zeros(N + 1)

for j in range(N + 1):
    # Interpolate to calculate the boundary layer thickness (y-component)
    d99[j] = np.interp(0.99*U[-1,j], U[:,j], y[:])
    
    # Calculate the displacement thickness
    integrand = 1 - U[:,j] / U[-1,j]
    d[j] = (dy / 2) * (integrand[0] + 2 * np.sum(integrand[1:-1]) + integrand[-1])

# Mean velocity Um and Vm
Um = np.zeros(N + 1)
Vm = np.zeros(N + 1)

for j in range(N + 1):
    Um[j] = np.mean(U[:,j])
    Vm[j] = np.mean(V[:,j])

print(
    f"\n\nThe solution has successfully converged.\n"
)

#%%
# Save the simulation as raw data
import os

print(
f"\nPlease wait while the data is being written to the output file.\n"
)

filename = "../data/simulation-data.csv"
directory = os.path.dirname(filename)  # Get the directory part of the file path

# Create the directory if it does not exist
if not os.path.exists(directory):
    os.makedirs(directory)

# Delete old file if it exists
if os.path.exists(filename):
    os.remove(filename)

# Open the file in append mode using built-in `open` instead of `os.open`
with open(filename, "w") as f:  # Use "w" mode to create a new file
    # Write the header
    f.write("x Um Vm d d99 tw\n")
    
    # Example loop to generate data
    for j in range(0, N + 1):
        f.write(f"{j} {Um[j]:.6f} {Vm[j]:.6f} {d[j]:.6f} {d99[j]:.6f} {T[0,j]:.6f}\n")
    f.write("#\n")
    
    for j in range(0, N + 1):
        f.write("x " + str(j) + "\n")
        f.write("y u v tau_xy\n")
        for i in range(0, M + 1):
            f.write(f"{i} {U[i,j]:.6f} {V[i,j]:.6f} {T[i,j]:.6f}\n")
        f.write("#\n")

# Open and read the file
# with open(filename, "r") as f:
#     print(f.read())
f.close()

print(
    f"The raw simulation data has been saved to the following file:\n"
    f"    {filename}\n"
)
#%%
# Plot the results

# import matplotlib

# matplotlib.rcParams['font.family'] = 'cmu serif'
# matplotlib.rcParams['mathtext.fontset'] = 'cm'
# matplotlib.rcParams['axes.unicode_minus'] = False
# plt = matplotlib.pyplot

# Plot velocity and shear stress @(k) [k as integer]
# k = 8       # s

# figure1, axis1 = plt.subplots(figsize = (4,3), dpi=300)

# lbl = r"$u(x=0,y)$"
# axis1.plot(U[:,0], y, label=lbl, color="grey", alpha=0.5, linestyle="--", linewidth="1")

# lbl = r"$u(x=" + str(k) + r",y)$"
# axis1.plot(U[:,k*10], y, label=lbl, color="black")
# axis1.plot([0, 0], [0, H], color="black")

# axis1.legend()
# axis1.set_ylim([0, 1])
# axis1.set_xlim([0, 1.1])

# axis1.set_ylabel("y")
# axis1.set_xlabel("u(x,y)")

# filename = "zpg_profile_zoom_x_" + str(k) + ".pdf"
# figure1.savefig(filename, bbox_inches='tight')


# figure2, axis2 = plt.subplots(figsize = (8,3), dpi=300)

# axis2.plot(U[:,0],y, color="black")
# axis2.plot([0,0], [0, H], color="black")

# axis2.plot(U[:,30] + 3, y, color="black")
# axis2.plot([3,3], [0, H], color="black")

# axis2.plot(U[:,60] + 6, y, color="black")
# axis2.plot([6,6], [0, H], color="black")

# axis2.plot(U[:,90] + 9, y, color="black")
# axis2.plot([9,9], [0, H], color="black")

# axis2.plot(x, d99, "--", linewidth=1)
# axis2.text(x[-1] + 0.25, d99[-1] - 0.025, r"$\mathbf{\delta}_{99}(x)$", color="#1f77b4")

# axis2.plot(x, d, "--", linewidth=1, color="#7f7f7f")
# axis2.text(x[-1] + 0.25, d[-1] - 0.025, r"$\mathbf{\delta}(x)$", color="#7f7f7f")

# axis2.set_ylim([0,2])
# axis2.set_xlim([-0.5, 11])
# axis2.set_xticks(np.arange(0, 11, 1))
# axis2.set_yticks(np.arange(0, 2.1, 0.5))

# axis2.set_ylabel("y")
# axis2.set_xlabel("x")
# axis2.set_title(r"Evolution of the velocity profile, $u(x,y)$")

# figure2.savefig("zpg_profile.pdf", bbox_inches='tight')


# figure3, axis3 = plt.subplots(figsize = (4,3), dpi=300)

# lbl = r"$\tau\:(x = 0,y)$"
# axis3.plot(T[:,0], y, label=lbl, color="grey", linestyle="--", linewidth=1)
# lbl = r"$\tau\:(x=" + str(k) + r",y)$"
# axis3.plot(T[:,k * 10], y, label=lbl, color="black")
# axis3.plot([0,0], [0,H], color="black")

# axis3.set_ylim([0,1])
# axis3.set_xlim([0,0.0055])
# axis3.set_xlabel(r"$\tau \:(x,y)$")
# axis3.set_ylabel("y")

# axis3.legend(loc="upper right")

# filename = "zpg_shear_zoom_x_" + str(k) + ".pdf"
# figure3.savefig(filename, bbox_inches='tight')


# figure4, axis4 = plt.subplots(figsize = (8,1.5), dpi = 300)

# axis4.plot(x, T[0,:], color="black")
# # axis4.text(x[-1] + 0.25, T[0, -1] - 0.00010, r"$\mathbf{\tau}_{w}(x)$", color="#1f77b4")

# axis4.set_ylim([0.002,0.005])
# axis4.set_xlim([-0,10])
# axis4.set_xticks(np.arange(0,11,1))

# axis4.set_ylabel(r"$\tau_{w}$")
# axis4.set_xlabel("x")
# axis4.set_title(r"Friction at the wall, $\tau_{w}(x)$")

# figure4.savefig("zpg_shear.pdf", bbox_inches='tight')