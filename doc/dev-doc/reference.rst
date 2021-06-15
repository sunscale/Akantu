.. _reference:

Reference
---------

Common
``````

.. doxygenfunction:: akantu::initialize(const std::string &input_file, int &argc, char **&argv)
.. doxygenfunction:: akantu::initialize(int &argc, char **&argv)

.. doxygentypedef:: akantu::UInt
.. doxygentypedef:: akantu::Int
.. doxygentypedef:: akantu::Real

.. doxygenenum:: akantu::ElementType
.. doxygenenum:: akantu::ModelType
.. doxygenenum:: akantu::AnalysisMethod
.. doxygenenum:: akantu::SolveConvergenceCriteria

.. doxygenclass:: akantu::ArrayBase
.. doxygenclass:: akantu::ArrayDataLayer
.. doxygenclass:: akantu::Array

.. doxygenclass:: akantu::ElementTypeMapArray

.. doxygenclass:: akantu::Vector
.. doxygenclass:: akantu::Matrix

Mesh
````
.. doxygenclass:: akantu::Mesh
.. doxygenclass:: akantu::FEEngine
.. doxygenclass:: akantu::Element
.. doxygenclass:: akantu::GroupManager
.. doxygenclass:: akantu::ElementGroup
.. doxygenclass:: akantu::NodeGroup

Models
``````

Common
......

.. doxygenclass:: akantu::BC::Dirichlet::FixedValue
.. doxygenclass:: akantu::BC::Dirichlet::FlagOnly
.. doxygenclass:: akantu::BC::Dirichlet::IncrementValue
.. doxygenclass:: akantu::BC::Neumann::FromHigherDim
.. doxygenclass:: akantu::BC::Neumann::FromSameDim
.. doxygenclass:: akantu::BoundaryCondition
.. doxygenclass:: akantu::BoundaryConditionFunctor
.. doxygenclass:: akantu::EventHandlerManager
.. doxygenclass:: akantu::Model
.. doxygenclass:: akantu::NonLocalManagerCallback

Solvers
.......

.. doxygenclass:: akantu::ModelSolver
.. doxygenclass:: akantu::DOFManager
.. doxygenclass:: akantu::NonLinearSolver
.. doxygenclass:: akantu::NonLinearSolverNewtonRaphson

Solid Mechanics Model
.....................

.. doxygenclass:: akantu::SolidMechanicsModel
.. doxygenclass:: akantu::SolidMechanicsModelOptions
.. doxygenclass:: akantu::MaterialSelector
.. doxygenclass:: akantu::MeshDataMaterialSelector
.. doxygenclass:: akantu::Material
.. doxygenclass:: akantu::InternalField

Solid Mechanics Model Cohesive
..............................

.. doxygenclass:: akantu::SolidMechanicsModelCohesive
.. doxygenclass:: akantu::FragmentManager

Heat Transfer Model
...................

.. doxygenclass:: akantu::HeatTransferModel

Phase Field Model
...................

.. doxygenclass:: akantu::PhaseFieldModel
.. doxygenclass:: akantu::PhaseField

Structural Mechanics Model
..........................

.. doxygenclass:: akantu::StructuralMaterial
.. doxygenclass:: akantu::StructuralMechanicsModel


Coupler Solid PhaseField
........................

.. doxygenclass:: akantu::CouplerSolidPhaseField


Synchronizers
`````````````
.. doxygenclass:: akantu::DataAccessor

Input/Output
````````````
.. doxygenclass:: akantu::Dumpable
.. doxygenclass:: akantu::DumperIOHelper
.. doxygenclass:: akantu::DumperParaview
.. doxygenclass:: akantu::DumperText
.. doxygenclass:: akantu::Field
.. doxygenclass:: akantu::Parser
.. doxygenclass:: akantu::ParserParameter
.. doxygenclass:: akantu::ParserSection
.. doxygenenum:: akantu::SectionType
