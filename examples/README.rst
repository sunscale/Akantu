========
Examples
========

This examples demonstrate some of the principal features of **Akantu**.  Each
sub-folder contain concentrate on one particular features.

Common features
---------------

*boundary_conditions*
  This shows how to impose boundary conditions by using the group of elements
  defining the boundaries. Boundaries are applied mainly by using functors, some
  are predefined in akantu asshown in the sub-example 'predefined_bc'. If this
  is not sufficient the user can define it's one functor, example
  'user_defined_bc'

*io*
  This examples show some advanced features of the dumpers, to dump extra
  fields than the default ones

*parallel*
  This example shows how to initialize a parallel computation

*python*
  This examples show some basic usage of the python interface

Solid Mechanics
---------------

*explicit*
  This examples shows a dynamic wave propagation in a beam. It uses a central
  difference explicit time integration scheme.

*implicit*
  This example shows how to use the implicit dynamic solver.

*static*
  This example shows how to do a static computation

*cohesive_element*
  This examples show some usage of the cohesive element in explicit or implicit
  and with extrinsic or intrinsic elements

*new_material*
  This example shows how to write a custom material outside of **Akantu** and
  how to use it

Heat transfer
-------------

*heat_transfer*
  This example shows how to use the HeatTransferModel

Structural mechanics
--------------------

*structural_mechanics*
  This example shows how to use the StructuralMechanicsModel
