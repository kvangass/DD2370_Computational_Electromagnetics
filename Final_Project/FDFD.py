#  We want to find the eigenvalue associated with the propagating mode
#  At certain frequencies there should only be a single propagating mode
"""
Code for the finite difference implementation
"""

import numpy as np
from scipy.sparse import diags, bmat, eye, save_npz
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import inv, eigs
import matplotlib.pyplot as plt
from scipy.io import savemat


def pos_to_index(nx, ny, Nx):
    m = (ny-1)*Nx + nx - 1  
    return m


def get_eps_grid(Nx2, Ny2, dx2, dy2, eps_sub, eps_core, eps_sup):
    """
    Generates the reference eps grid

    n_sup = 1
    n_core = 1.9
    n_sub = 1.52

    eps_sup = 1
    eps_core = n_core**2
    eps_sub = n_sub**2
    """
    eps_grid2 = np.ones((Ny2, Nx2))*eps_sup
    # implement rib
    rib_width = 2e-6
    rib_width = round(rib_width / dx2 / 2)
    rib_height = round(0.6e-6 / dy2 / 2)
    rib_h_pos = round((6e-6 - (2 + 0.25 + 0.3) * 1e-6) / dy2)
    eps_grid2[rib_h_pos - rib_height: rib_h_pos + rib_height,
              Nx - rib_width:Nx + rib_width] = eps_core
    # implement lower core
    l_core_height = 0.25e-6
    l_core_height = round(l_core_height / dy2 / 2)

    eps_grid2[rib_h_pos + rib_height: rib_h_pos + rib_height + l_core_height,
              :] = eps_core
    # implement sub
    sub_height = rib_h_pos + rib_height + l_core_height
    eps_grid2[sub_height:, :] = eps_sub
    return eps_grid2


def get_eps_circ(Nx2, Ny2, dx2, dy2, eps_core, eps_clad, r_core):
    """
    Generates the epsilons for a circular waveguide
    """
    eps_grid2 = np.ones((Ny2, Nx2)) * eps_clad
    x_center = int(Nx2 / 2) * dx2
    y_center = int(Ny2 / 2) * dy2

    mask = np.zeros(eps_grid2.shape) == 1
    for x_val in range(Nx2):
        for y_val in range(Ny2):
            x_pos = dx2 * x_val
            y_pos = dx2 * y_val

            r_squared = (x_pos - x_center) ** 2 + (y_pos - y_center) ** 2

            in_circle = r_squared < r_core ** 2
            if in_circle:
                mask[y_val, x_val] = True
    eps_grid2[mask] = eps_core
    return eps_grid2


def diff_matrices(Nx, Ny, dx, dy, k0):
    """
    Generates the differential matrices for D_e_x , D_e_y, D_h_x, D_h_y
    m = (ny - 1)*Nx + nx
    Using column order grid
    """
    dx = dx*k0  
    dy = dy*k0
    N = Nx*Ny
    diagonals_x = [-1/dx*np.ones(N), 1/dx*np.ones(N)]
    offsets_x = [0, 1]
    D_e_x = diags(diagonals_x, offsets_x)
    D_e_x = D_e_x.tolil()
    #Implement Dirchlet BC
    for ny in range(Ny):
        # For points on the x boundary set the
        # relation to the next point to zero
        idx1 = pos_to_index(Nx, ny, Nx)  # This node
        idx2 = pos_to_index(Nx+1, ny, Nx)  # Set This relation
        D_e_x[idx1, idx2] = 0  # Most cases are already 0
    D_e_x = D_e_x.tocsc()
    diagonals_y = [-1/dy*np.ones(N), 1/dy*np.ones(N)]
    offsets_y = [0, Nx]
    D_e_y = diags(diagonals_y, offsets_y)
    D_e_y = D_e_y.tocsc()
    D_h_x = -D_e_x.getH()  # May need to change the type of Dex to improve speed
    D_h_y = -D_e_y.getH()
    return D_e_x, D_e_y, D_h_x, D_h_y


## CONSTANTS
order_str = 'C'  # Gives the unravel order of the grid
eps0 = 8.8541878e-12
mu0 = 4e-7 * np.pi
c0 = 299792458

#freq_0 = 1.5e9 #1.157e9
#k0 = 2*np.pi*freq_0/c0
lambda0 = 3e-6 #0.9e-6
k0 = 2*np.pi/lambda0

n_clad = 1.0
n_core = 3.0 #1.9

eps_clad = n_clad**2
eps_core = n_core**2

# Cell size
#h = 0.0333e-6
h = 0.02e-6
# Cell sizes smaller than h=0.033e-6 take way too long.
# Waveguide Dimensions
# Cutoff at 1.157 GHz
Lx = 6e-6  # um
Ly = 6e-6  # um

r_core = 1e-6

# Number of cells on x and y axis
# This is snap to grid for the waveguide
Nx = int(np.round(Lx / h))
Ny = int(np.round(Ly / h))

print(f'N: {Nx*Ny}')

Nx2 = Nx*2
Ny2 = Ny*2

dx = Lx / Nx
dy = Ly / Ny

dx2 = Lx / Nx2
dy2 = Ly / Ny2

# Number of eigenmodes to solve
N_eig = 4

# Creates the double grid shape representation
eps_grid2 = np.ones((Ny2, Nx2))*eps_clad
x_center = int(Nx2/2)
y_center = int(Ny/2)


#eps_grid2[core_mask] = eps_core
eps_grid2 = get_eps_circ(Nx2, Ny2, dx2, dy2, eps_core, eps_clad, r_core)
mu_grid2 = np.ones(eps_grid2.shape)  # this may not be necessary

# Numpy flattens row wise

# ERxx = ER2(2:2:Nx2,1:2:Ny2)
eps_xx2 = eps_grid2[1:Nx2:2, 0:Ny2:2]
eps_xx = diags(eps_xx2.flatten(order=order_str))  # order='F' is fortran

# ERyy = ER2(1:2:Nx2, 2:2:Ny2)
eps_yy2 = eps_grid2[0:Nx2:2, 1:Ny2:2]
eps_yy = diags(eps_yy2.flatten(order=order_str))

# ERzz = ER2(1:2:Nx2, 1:2:Ny2)
eps_zz2 = eps_grid2[0:Nx2:2, 0:Ny2:2]
eps_zz = diags(eps_zz2.flatten(order=order_str))
inv_eps_zz = inv(eps_xx.tocsc())

# URxx = UR2(1:2:Nx2, 2:2:Ny2)
mu_xx2 = mu_grid2[0:Nx2:2, 1:Ny2:2]
mu_xx = diags(mu_xx2.flatten(order=order_str))

# URyy = UR2(2:2:Nx2, 1:2:Ny2)
mu_yy2 = mu_grid2[1:Nx2:2, 0:Ny2:2]
mu_yy = diags(mu_yy2.flatten(order=order_str))

# URzz = UR2(2:2:Nx2, 2:2:Ny2)
mu_zz2 = mu_grid2[1:Nx2:2, 1:Ny2:2]
mu_zz = diags(mu_zz2.flatten(order=order_str))
inv_mu_zz = inv(mu_zz.tocsc())

##########################
D_e_x, D_e_y, D_h_x, D_h_y = diff_matrices(Nx, Ny, dx, dy, k0)

p_xx = D_e_x.dot(inv_eps_zz).dot(D_h_y)
p_xy = -(D_e_x.dot(inv_eps_zz).dot(D_h_x) + mu_yy)
p_yx = D_e_y.dot(inv_eps_zz).dot(D_h_y) + mu_xx
p_yy = -D_e_y.dot(inv_eps_zz).dot(D_h_x)
P = bmat([[p_xx, p_xy], [p_yx, p_yy]])

q_xx = D_h_x.dot(inv_mu_zz).dot(D_e_y)
q_xy = -(D_h_x.dot(inv_mu_zz).dot(D_e_x) + eps_yy)
q_yx = D_h_y.dot(inv_mu_zz).dot(D_e_y) + eps_xx
q_yy = -D_h_y.dot(inv_mu_zz).dot(D_e_x)
Q = bmat([[q_xx, q_xy], [q_yx, q_yy]])

Omega = P.dot(Q)

ex_conductor_grid = np.zeros((Ny, Nx))
ey_conductor_grid = np.zeros((Ny, Nx))

# Boundary 1 and 3: Ex = 0
# Top and bottom of rectangular wg
ex_conductor_grid[0, :] = 1
ex_conductor_grid[Ny-1, :] = 1
#ex_conductor_grid[:, 0] = 1
#ex_conductor_grid[:, Nx-1] = 1
# Boundary 2 and 4: Ey = 0
# Left and Right of rectangular wg
#ey_conductor_grid[0, :] = 1
#ey_conductor_grid[Ny-1, :] = 1
ey_conductor_grid[:, 0] = 1
ey_conductor_grid[:, Nx-1] = 1

# Generate Force matrix: F
# The Force matrix shows which rows of omega to zero
# Flatten and combine grids
flat_ex_conductor = ex_conductor_grid.flatten(order_str)
flat_ey_conductor = ey_conductor_grid.flatten(order_str)

F_vec = np.zeros((Nx*Ny*2))
F_vec[0:Nx*Ny] = flat_ex_conductor
F_vec[Nx*Ny::] = flat_ey_conductor

F = diags(F_vec, 0)

# sparse eye
I_matrix = eye(Nx*Ny*2)

Omega_d_prime = F + (I_matrix - F).dot(Omega)
# Solve the eigenvalue problem
# v_omega are the eigenvalues, gamma_tilde**2
# w_omega are the eigenvectors
# Need to relate the gamma_tilde**2 to the actual eigenvalues
v_omega, w_omega = eigs(Omega_d_prime, k=N_eig, sigma=-n_core**2)

# Select out the components of the eigenvector
#kz_solved = np.imag(-np.sqrt(v_omega)*k0)  # Not sure about this algebra

if order_str == 'F':
    fields = np.reshape(w_omega, (Ny, Nx, 2, N_eig), order=order_str)
elif order_str == 'C':
    w_omega = w_omega.transpose()
    fields = np.reshape(w_omega, (N_eig, 2, Nx, Ny), order=order_str)

n_eff = np.imag(np.sqrt(v_omega))
print('N_eff:')
print(n_eff)

# Figures
plt.figure(50)
plt.title('eps_grid')
plt.imshow(eps_grid2)


for eig_index in range(N_eig):
    plt.figure(eig_index)
    plt.title(f'$|E_x|$ Mode {eig_index}')
    if order_str == 'F':
        plt.imshow(np.abs(fields[:, :, 0, eig_index]))
    elif order_str == 'C':
        plt.imshow(np.abs(fields[eig_index, 0, :, :]))
    plt.clim(0, 0.1)
    plt.colorbar()
    plt.savefig(f'./Final_report/Figures/e_x_mode_{eig_index}.png')

for eig_index in range(N_eig):
    plt.figure(eig_index+N_eig)
    plt.title(f'|$E_y|$ Mode {eig_index}')
    if order_str == 'F':
        plt.imshow(np.abs(fields[:, :, 1, eig_index]))
    elif order_str == 'C':
        plt.imshow(np.abs(fields[eig_index, 1, :, :]))
    plt.clim(0, 0.1)
    plt.colorbar()
    plt.savefig(f'./Final_report/Figures/e_y_mode{eig_index}.png')

plt.show()

