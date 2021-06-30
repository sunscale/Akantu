Version 4.0
-------------

features
    * transferred the CI from jenkins to Gitlab CI/CD

api
    * breaking api changes to match the stl containers

      * ``clear`` does not set to ``0`` anymore but empties containers
      * ``empty`` does not empty containers but tells if the container is empty
      * ``zero`` replace the old empty and set containers to ``0``

Version 3.2
-------------

features
    Activating PETSc solver back with the new solver interface

api
    deprecating old C++ 03 code

Version 3.0
-------------

c++14
    switch from C++ standard ``2003`` to ``2014``
    Example of changes implied by this::

      for (UInt g = _not_ghost; g <= _ghost; ++g) {
        GhostType gt = (GhostType)g;
        Mesh::type_iterator it = this->mesh.firstType(spatial_dimension, gt);
        Mesh::type_iterator end = this->mesh.lastType(spatial_dimension, gt);
        for (; it != end; ++it) {
          ElementType & type = *it;
          ...
        }
      }

    becomes::

      for (auto ghost_type : ghost_types) {
        for (auto type : mesh.elementTypes(spatial_dimension,
                                           ghost_type)) {
          ...
        }
      }

features
   * Parallel cohesive elements
   * Models using new interface for solvers

      * Same configuration for all models
      * Solver can be configured in input file
      * PETSc interface temporary inactive

   * Periodic boundary condition temporary inactive
   * Element groups created by default for ``“physical_names”``

api
   * Named arguments for functions (e.g. ``model.initFull(_analysis_method = _static)``)
   * Only one function to solve a step :cpp:func:`model.solveStep() <akantu::ModelSolver::solveStep>`
   * Simplification of the parallel simulation with the :cpp:func:`mesh.distribute() <akantu::Mesh::distribute>` function
