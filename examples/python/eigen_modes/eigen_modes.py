#!/usr/bin/env python

import subprocess
import argparse
import akantu as aka
import numpy as np
from scipy.sparse.linalg import eigsh
from scipy.sparse import csr_matrix
try:
    import matplotlib.pyplot as plt
    from image_saver import ImageSaver
    has_matplotlib = True
except ImportError:
    has_matplotlib = False

# -----------------------------------------------------------------------------
# parser
# -----------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Eigen mode exo')
parser.add_argument('-m', '--mode_number', type=int,
                    help='precise the mode to study', default=2)

parser.add_argument('-wL', '--wave_width', type=float,
                    help='precise the width of the wave for '
                    'the initial displacement', default=5)

parser.add_argument('-L', '--Lbar', type=float,
                    help='precise the length of the bar', default=10)

parser.add_argument('-t', '--time_step', type=float,
                    help='precise the timestep',
                    default=None)

parser.add_argument('-n', '--max_steps', type=int,
                    help='precise the number of timesteps',
                    default=500)

parser.add_argument('-mh', '--mesh_h', type=float,
                    help='characteristic mesh size',
                    default=.2)

parser.add_argument('-p', '--plot', action='store_true',
                    help='plot the results')

args = parser.parse_args()
mode = args.mode_number
wave_width = args.wave_width
time_step = args.time_step
max_steps = args.max_steps
mesh_h = args.mesh_h
Lbar = args.Lbar
plot = args.plot

# -----------------------------------------------------------------------------
# Mesh Generation
# -----------------------------------------------------------------------------
geo_content = """
// Mesh size
h  = {0};
""".format(mesh_h)

geo_content += """
h1 = h;
h2 = h;

// Dimensions of the bar
Lx = 10;
Ly = 1;

// ------------------------------------------
// Geometry
// ------------------------------------------

Point(101) = { 0.0, -Ly/2, 0.0, h1};
Point(102) = { Lx,  -Ly/2, 0.0, h2};

Point(103) = { Lx,  0., 0.0,  h2};
Point(104) = { Lx,  Ly/2., 0.0,  h2};

Point(105) = { 0.0, Ly/2., 0.0,  h1};
Point(106) = { 0.0, 0., 0.0,  h1};

Line(101) = {101,102};
Line(102) = {102,103};
Line(103) = {103,104};
Line(104) = {104,105};
Line(105) = {105,106};
Line(106) = {106,101};
Line(107) = {106,103};


Line Loop(108) = {101, 102, -107, 106};
Plane Surface(109) = {108};
Line Loop(110) = {103, 104, 105, 107};
Plane Surface(111) = {110};
Physical Surface(112) = {109, 111};

Transfinite Surface "*";
Recombine Surface "*";
Physical Surface(113) = {111, 109};

Physical Line("XBlocked") = {103, 102};
Physical Line("ImposedVelocity") = {105, 106};
Physical Line("YBlocked") = {104, 101};
"""

mesh_file = 'bar'
with open(mesh_file + '.geo', 'w') as f:
    f.write(geo_content)

subprocess.call(['gmsh', '-2', mesh_file + '.geo'])
mesh_file = mesh_file + '.msh'


# -----------------------------------------------------------------------------
# Initialization
# -----------------------------------------------------------------------------
spatial_dimension = 2
aka.parseInput('material.dat')

mesh = aka.Mesh(spatial_dimension)
mesh.read(mesh_file)

model = aka.SolidMechanicsModel(mesh)
model.initFull(aka._implicit_dynamic)

model.setBaseName("waves-{0}".format(mode))
model.addDumpFieldVector("displacement")
model.addDumpFieldVector("acceleration")
model.addDumpFieldVector("velocity")
model.addDumpField("blocked_dofs")

# -----------------------------------------------------------------------------
# Boundary conditions
# -----------------------------------------------------------------------------
internal_force = model.getInternalForce()
displacement = model.getDisplacement()
acceleration = model.getAcceleration()
velocity = model.getVelocity()

blocked_dofs = model.getBlockedDOFs()
nbNodes = mesh.getNbNodes()
position = mesh.getNodes()

model.applyBC(aka.FixedValue(0.0, aka._x), "XBlocked")
model.applyBC(aka.FixedValue(0.0, aka._y), "YBlocked")

# ------------------------------------------------------------------------
# timestep value computation
# ------------------------------------------------------------------------

time_factor = 0.8
stable_time_step = model.getStableTimeStep() * time_factor

if time_step:
    print("Required Time Step = {0}".format(time_step))
    if stable_time_step * time_factor < time_step:
        print("Stable Time Step = {0}".format(stable_time_step))
        raise RuntimeError("required time_step too large")
    print("Required Time Step = {0}".format(time_step))
else:
    print("Stable Time Step = {0}".format(stable_time_step))
    time_step = stable_time_step * time_factor

model.setTimeStep(time_step)


# ------------------------------------------------------------------------
# compute the eigen modes
# ------------------------------------------------------------------------
model.assembleStiffnessMatrix()
model.assembleMass()
stiff = model.getDOFManager().getMatrix('K')
stiff = aka.AkantuSparseMatrix(stiff).toarray()
mass = model.getDOFManager().getMatrix('M')
mass = aka.AkantuSparseMatrix(mass).toarray()

# select the non blocked DOFs by index in the mask
mask = np.equal(blocked_dofs.flatten(), False)

Mass_star = mass[mask, :]
Mass_star = csr_matrix(Mass_star[:, mask].copy())

K_star = stiff[mask, :]
K_star = csr_matrix(K_star[:, mask].copy())

print('getting the eigen values')
vals, vects = eigsh(K_star, M=Mass_star, which='SM', k=20)

# -----------------------------------------------------------------------------
# import the initial conditions in displacement
# -----------------------------------------------------------------------------
displacement.reshape(nbNodes*2)[mask] = vects[:, mode]
with open('modes.txt', 'a') as f:
    f.write('{0} {1}\n'.format(mode, vals[mode]))

model.dump()

# -----------------------------------------------------------------------------
# prepare the storage of the dynamical evolution
# -----------------------------------------------------------------------------
e_p = np.zeros(max_steps + 1)
e_k = np.zeros(max_steps + 1)
e_t = np.zeros(max_steps + 1)
time = np.zeros(max_steps + 1)
norm = np.zeros(max_steps + 1)

epot = model.getEnergy('potential')
ekin = model.getEnergy('kinetic')

e_p[0] = epot
e_k[0] = ekin
e_t[0] = epot + ekin
time[0] = 0

if has_matplotlib:
    disp_sav = ImageSaver(mesh, displacement, 0, Lbar)
    velo_sav = ImageSaver(mesh, velocity, 0, Lbar)


# -----------------------------------------------------------------------------
# loop for evolution of motion dynamics
# -----------------------------------------------------------------------------
for step in range(1, max_steps + 1):
    model.solveStep()
    # outputs
    epot = model.getEnergy('potential')
    ekin = model.getEnergy('kinetic')

    print(step, '/', max_steps, epot, ekin, epot + ekin)

    e_p[step] = epot
    e_k[step] = ekin
    e_t[step] = epot + ekin
    time[step] = (step + 1) * time_step

    if has_matplotlib:
        disp_sav.storeStep()
        velo_sav.storeStep()

    if step % 10 == 0:
        model.dump()


if plot and has_matplotlib:
    # --------------------------------------------------------------------------
    # plot figures for global evolution
    # --------------------------------------------------------------------------
    # energy norms
    plt.figure(1)
    plt.plot(time, e_t, 'r', time, e_p, 'b', time, e_k, 'g')

    # space-time diagram for diplacements
    plt.figure(2)
    plt.imshow(disp_sav.getImage(), extent=(0, Lbar, max_steps * time_step, 0))
    plt.xlabel("Space ")
    plt.ylabel("Time ")

    # space-time diagram for velocities
    plt.figure(3)
    plt.imshow(velo_sav.getImage(), extent=(0, Lbar, max_steps * time_step, 0))
    plt.xlabel("Velocity")
    plt.ylabel("Time")
    plt.show()
