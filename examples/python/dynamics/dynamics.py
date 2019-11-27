#!/usr/bin/env python3
import numpy as np
import akantu as aka


# -----------------------------------------------------------------------------
class MyFixedValue(aka.FixedValue):
    def __init__(self, value, axis):
        super().__init__(value, axis)
        self.value = value
        self.axis = int(axis)

    def __call__(self, node, flags, disp, coord):
        # sets the displacement to the desired value in the desired axis
        disp[self.axis] = self.value
        # sets the blocked dofs vector to true in the desired axis
        flags[self.axis] = True


# -----------------------------------------------------------------------------
def main():
    spatial_dimension = 2
    mesh_file = 'bar.msh'
    max_steps = 250
    time_step = 1e-3

    aka.parseInput('material.dat')

    # -------------------------------------------------------------------------
    # Initialization
    # -------------------------------------------------------------------------
    mesh = aka.Mesh(spatial_dimension)
    mesh.read(mesh_file)

    model = aka.SolidMechanicsModel(mesh)

    model.initFull(_analysis_method=aka._explicit_lumped_mass)
    # model.initFull(_analysis_method=aka._implicit_dynamic)

    model.setBaseName("waves")
    model.addDumpFieldVector("displacement")
    model.addDumpFieldVector("acceleration")
    model.addDumpFieldVector("velocity")
    model.addDumpFieldVector("internal_force")
    model.addDumpFieldVector("external_force")
    model.addDumpField("strain")
    model.addDumpField("stress")
    model.addDumpField("blocked_dofs")

    # -------------------------------------------------------------------------
    # boundary conditions
    # -------------------------------------------------------------------------
    model.applyBC(MyFixedValue(0, aka._x), "XBlocked")
    model.applyBC(MyFixedValue(0, aka._y), "YBlocked")

    # -------------------------------------------------------------------------
    # initial conditions
    # -------------------------------------------------------------------------
    displacement = model.getDisplacement()
    nb_nodes = mesh.getNbNodes()
    position = mesh.getNodes()

    pulse_width = 1
    A = 0.01
    for i in range(0, nb_nodes):
        # Sinus * Gaussian
        x = position[i, 0] - 5.
        L = pulse_width
        k = 0.1 * 2 * np.pi * 3 / L
        displacement[i, 0] = A * \
            np.sin(k * x) * np.exp(-(k * x) * (k * x) / (L * L))
        displacement[i, 1] = 0

    # -------------------------------------------------------------------------
    # timestep value computation
    # -------------------------------------------------------------------------
    time_factor = 0.8
    stable_time_step = model.getStableTimeStep() * time_factor

    print("Stable Time Step = {0}".format(stable_time_step))
    print("Required Time Step = {0}".format(time_step))

    time_step = stable_time_step * time_factor

    model.setTimeStep(time_step)

    # -------------------------------------------------------------------------
    # loop for evolution of motion dynamics
    # -------------------------------------------------------------------------
    print("step,step * time_step,epot,ekin,epot + ekin")
    for step in range(0, max_steps + 1):

        model.solveStep()

        if step % 10 == 0:
            model.dump()

        epot = model.getEnergy('potential')
        ekin = model.getEnergy('kinetic')

        # output energy calculation to screen
        print("{0},{1},{2},{3},{4}".format(step, step * time_step,
                                           epot, ekin,
                                           (epot + ekin)))

    return


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
