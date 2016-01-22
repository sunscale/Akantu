import akantu
import numpy as np

################################################################
## Dirichlet Boudary condition functor: fix the displacement
################################################################
class FixedValue:

    def __init__(self,value,axis):
        self.value = value
        if axis == 'x': axis = 0
        if axis == 'y': axis = 1
        self.axis = axis
        
    def operator(self,node,flags,disp,coord):
        # sets the displacement to the desired value in the desired axis
        disp[self.axis] = self.value
        # sets the blocked dofs vector to true in the desired axis
        flags[self.axis] = True

################################################################
## Neumann Boudary condition functor: impose a traction
################################################################
class FromTraction:

    def __init__(self,traction):
        self.traction = traction

    def operator(self,quad_point,force,coord,normals):
        # sets the force to the desired value in the desired axis
        force[:] = self.traction

################################################################

def solve(material_file,mesh_file,traction):
    akantu.initialize(material_file)
    spatial_dimension = 2
    
    ################################################################
    ##Initialization                                                          
    ################################################################
    mesh = akantu.Mesh(spatial_dimension)
    mesh.read(mesh_file)
    akantu.MeshUtils.purifyMesh(mesh)
    mesh.createGroupsFromStringMeshData("physical_names")

    model = akantu.SolidMechanicsModel(mesh)
    model.initFull(akantu.SolidMechanicsModelOptions(akantu._static))
    model.assembleStiffnessMatrix()
    model.updateResidual()
    
    model.setBaseName("plate")
    model.addDumpFieldVector("displacement")
    model.addDumpFieldVector("force")
    model.addDumpField("boundary")
    model.addDumpField("strain")
    model.addDumpField("stress")
    model.addDumpField("blocked_dofs")
    
    ################################################################
    ##Boundary conditions                                                      
    ################################################################
    displacement = model.getDisplacement()
    blocked_dofs = model.getBlockedDOFs()
    
    model.applyDirichletBC(FixedValue(0.0, 'x'), "XBlocked")
    model.applyDirichletBC(FixedValue(0.0, 'y'), "YBlocked")

    trac = np.zeros(spatial_dimension)
    trac[1] = traction

    print "Solve for traction ",traction

    model.getForce()[:] = 0
    model.applyNeumannBC(FromTraction(trac), "Traction")

    model.solveStaticDisplacement(1e-10,2);

    model.dump()
    akantu.finalize()

################################################################
# main 
################################################################

def main():

    import os
    mesh_file = 'plate.msh'
    #if mesh was not created the calls gmsh to generate it
    if not os.path.isfile(mesh_file):
        import subprocess
        subprocess.call('gmsh -2 plate.geo {0}'.format(mesh_file),shell=True)
    
    material_file = 'material.dat'
    spatial_dimension = 2

    traction = 1.
    solve(material_file,mesh_file,traction)

################################################################
if __name__ == "__main__":
    main() 

