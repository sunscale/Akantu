#!/usr/bin/python
################################################################
import os
import subprocess
import numpy as np
import akantu
################################################################

class FixedValue:

    def __init__(self, value, axis):
        self.value = value
        if axis == 'x':
            axis = 0
        if axis == 'y':
            axis = 1
        self.axis = axis

    def operator(self, node, flags, disp, coord):
        # sets the displacement to the desired value in the desired axis
        disp[self.axis] = self.value
        # sets the blocked dofs vector to true in the desired axis
        flags[self.axis] = True

################################################################


class LocalElastic:

    def __init__(self):

        ## young modulus
        self.E = 1
        ## Poisson coefficient
        self.nu = 0.3
        ## density
        self.rho = 1
        ## First Lame coefficient
        self._lambda = self.nu * self.E / ((1 + self.nu) * (1 - 2*self.nu))
        ## Second Lame coefficient (shear modulus)
        self.mu = self.E / (2 * (1 + self.nu));

    ## declares all the internals
    def registerInternals(self):
        return []

    ## declares all the parameters that could be parsed
    def registerParam(self):
        return []


    ## declares all the parameters that are needed
    def getPushWaveSpeed(self):
        return np.sqrt((self._lambda + 2*self.mu)/self.rho);

    ## constitutive law for a given quadrature point
    def computeStress(self,grad_u,sigma,internals):
        lbda = 1.
        mu = 1.

        trace = grad_u.trace()
        sigma[:, :] = lbda*trace*np.eye(2) + mu * (grad_u + grad_u.T)


################################################################
def main():

    spatial_dimension = 2
    Lbar = 10.
    akantu.initialize('material.dat')

    mesh_file = 'bar.msh'
    max_steps = 250
    time_step = 1e-3

    #if mesh was not created the calls gmsh to generate it
    if not os.path.isfile(mesh_file):
        ret = subprocess.call('gmsh -2 bar.geo bar.msh', shell=True)
        if ret != 0:
            raise Exception('execution of GMSH failed: do you have it installed ?')


    ################################################################
    ## Initialization
    ################################################################
    mesh = akantu.Mesh(spatial_dimension)
    mesh.read(mesh_file)

    mesh.createGroupsFromStringMeshData("physical_names")
    model = akantu.SolidMechanicsModel(mesh)

    model.initFull(akantu.SolidMechanicsModelOptions(akantu._explicit_lumped_mass, True))
    mat = LocalElastic()
    model.registerNewPythonMaterial(mat, "local_elastic")
    model.initMaterials()


    model.setBaseName("waves")
    model.addDumpFieldVector("displacement")
    model.addDumpFieldVector("acceleration")
    model.addDumpFieldVector("velocity")
    model.addDumpField("blocked_dofs")

    ################################################################
    ## Boundary conditions
    ################################################################
    residual = model.getResidual()
    mass = model.getMass()

    displacement = model.getDisplacement()
    acceleration = model.getAcceleration()
    velocity = model.getVelocity()

    blocked_dofs = model.getBlockedDOFs()

    ################################################################
    ## boundary conditions
    ################################################################

    model.applyDirichletBC(FixedValue(0, 'x'), "XBlocked")
    model.applyDirichletBC(FixedValue(0, 'y'), "YBlocked")

    ################################################################
    ## initial conditions
    ################################################################

    nb_nodes = mesh.getNbNodes()
    position = mesh.getNodes()

    pulse_width = 1
    A = 0.01
    for i in range(0, nb_nodes):
        # Sinus * Gaussian
        x = position[i, 0] - 5.
        L = pulse_width
        k = 0.1 * 2 * np.pi * 3 / L
        displacement[i, 0] = A * np.sin(k * x) * np.exp(-(k * x) * (k * x) / (L * L))

    ################################################################
    ## timestep value computation
    ################################################################
    time_factor = 0.8
    stable_time_step = model.getStableTimeStep() * time_factor

    print "Stable Time Step = {0}".format(stable_time_step)
    print "Required Time Step = {0}".format(time_step)

    time_step = stable_time_step * time_factor

    model.setTimeStep(time_step)

    ################################################################
    ## loop for evolution of motion dynamics
    ################################################################
    model.updateResidual()

    epot = model.getEnergy('potential')
    ekin = model.getEnergy('kinetic')

    print "step,step * time_step,epot,ekin,epot + ekin"
    for step in range(0, max_steps+1):

        model.dump()
        ## output energy calculation to screen
        print "{0},{1},{2},{3},{4}".format(step, step * time_step,
                                           epot, ekin,
                                           (epot + ekin))

        model.solveStep()

    akantu.finalize()
    return

################################################################
if __name__ == "__main__":
    main()
