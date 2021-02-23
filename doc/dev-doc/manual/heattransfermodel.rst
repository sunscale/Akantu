Heat Transfer Model
===================

The heat transfer model is a specific implementation of the :cpp:class:`Model
<akantu::Model>` interface dedicated to handle the dynamic heat equation.

Theory
------

The strong form of the dynamic heat equation
can be expressed as

.. math::
  \rho c_v \dot{T} + \nabla \cdot \vec{\kappa} \nabla T = b

with :math:`T` the scalar temperature field, :math:`c_v` the specific heat capacity, :math:`\rho`
the mass density, :math:`\mat{\kappa}` the conductivity tensor, and :math:`b` the heat
generation per unit of volume. The discretized weak form with a finite number of
elements is

.. math::
  \forall i \quad
  \sum_j \left( \int_\Omega \rho c_v N_j N_i  d\Omega \right) \dot{T}_j
  - \sum_j \left( \int_\Omega \vec{\kappa} \nabla N_j \nabla N_i d\Omega \right) T_j =
  - \int_{\Gamma}  N_i \vec{q} \cdot \vec{n} d\Gamma + \int_\Omega b N_i d\Omega

with :math:`i` and :math:`j` the node indices, :math:`\vec{n}` the normal field to the surface
:math:`\Gamma = \partial \Omega`.
To simplify, we can define the capacity and the conductivity matrices as

.. math::
  C_{ij} = \int_\Omega \rho c_v N_j N_i  d\Omega \qquad \textrm{and} \qquad
  K_{ij} = - \int_\Omega \vec{\kappa} \nabla N_j \nabla N_i d\Omega

and the system to solve can be written

.. math::
  \mat{C} \cdot \vec{\dot{T}} = \vec{Q}^{\text{ext}} -\mat{K} \cdot \vec{T}~,

with :math:`\vec{Q}^{\text{ext}}` the consistent heat generated.

Using the Heat Transfer Model
-----------------------------

A material file name has to be provided during initialization.
Currently, the :cpp:class:`HeatTransferModel <akantu::HeatTransferModel>` object uses dynamic analysis
with an explicit time integration scheme.  It can simply be created
like this

.. code-block:: c++

   HeatTransferModel model(mesh, spatial_dimension);

while an existing mesh has been used (see \ref{sect:common:mesh}).
Then the model object can be initialized with:

.. code-block:: c++

   model.initFull()

This function will load the material properties, and allocate / initialize the nodes and element :cpp:class:`Arrays <akantu::Array>`
More precisely, the heat transfer model contains 4 :cpp:class:`Arrays <akantu::Array>`:

- **temperature** contains the nodal temperature :math:`T` (zero by default after the initialization).

- **temperature_rate** contains the variations of temperature :math:`\dot{T}` (zero by default after the initialization).

- **blocked_dofs** contains a Boolean value for each degree of freedom specifying whether the degree is blocked or not. A Dirichlet boundary condition (:math:`T_d`) can be prescribed by setting the **blocked_dofs** value of a degree of freedom to ``true``. The **temperature** and the **temperature_rate** are computed for all degrees of freedom where the **blocked_dofs** value is set to ``false``. For the remaining degrees of freedom, the imposed values (zero by default after initialization) are kept.

- **external_heat_rate** contains the external heat generations. :math:`\vec{Q^{ext}}` on the nodes.

- **internal_heat_rate** contains the internal heat generations. :math:`\vec{Q^{int}} = -\mat{K} \cdot \vec{T}` on the nodes.

Only a single material can be specified on the domain. A material text file (*e.g.* material.dat) provides the material properties as follows:

.. code-block::

  model heat_transfer_model [
    capacity = %\emph{XXX}%
    density = %\emph{XXX}%
    conductivity = [%\emph{XXX}% ... %\emph{XXX}%]
  ]

where the ``capacity`` and ``density`` are scalars, and the ``conductivity`` is specified as a :math:`3\times 3` tensor.

Explicit Dynamic
----------------

The explicit  time integration scheme in ``Akantu``  uses a lumped capacity
matrix :math:`\mat{C}` (reducing the computational  cost, see Chapter :ref:`sect-smm`).
This matrix is assembled by distributing the capacity of each element onto its nodes. Therefore, the resulting :math:`\mat{C}` is a diagonal matrix stored in the ``capacity`` :cpp:class:`Array <akantu::Array>` of the model.


.. code-block:: c++

   model.assembleCapacityLumped();

.. note::
   Currently, only the explicit time integration with lumped capacity
   matrix is implemented within ``Akantu``.

The explicit integration scheme is *Forward Euler* :cite:`curnier92a`.

- Predictor: :math:`\vec{T}_{n+1} = \vec{T}_{n} + \Delta t \dot{\vec{T}}_{n}`

- Update residual: :math:`\vec{R}_{n+1} = \left( \vec{Q^{ext}_{n+1}} - \vec{K}\vec{T}_{n+1} \right)`

- Corrector : :math:`\dot{\vec{T}}_{n+1} = \mat{C}^{-1} \vec{R}_{n+1}`

The explicit integration scheme is conditionally stable. The time step has to be
smaller than the stable time step, and it can be obtained in ``Akantu`` as
follows:

.. code-block:: c++

   time_step = model.getStableTimeStep();

The stable time step is defined as:

.. math::
  \Delta t_{\st{crit}} = 2 \Delta x^2 \frac{\rho c_v}{\mid\mid \mat{\kappa} \mid\mid^\infty}
  :label: eqn:htm:explicit:stabletime

where :math:`\Delta x` is the characteristic length (*e.g* the in-radius in the
case of linear triangle element), :math:`\rho` is the density,
:math:`\mat{\kappa}` is the conductivity tensor, and :math:`c_v` is the specific
heat capacity. It is necessary to impose a time step which is smaller than the
stable time step, for instance, by multiplying the stable time step by a safety
factor smaller than one.

.. code-block:: c++

   const Real safety_time_factor = 0.1;
   Real applied_time_step = time_step * safety_time_factor;
   model.setTimeStep(applied_time_step);


The following loop allows, for each time step, to update the ``temperature``,
``residual`` and ``temperature_rate`` fields following the previously described
integration scheme.

.. code-block:: c++

   for (UInt s = 1; (s-1)*applied_time_step < total_time; ++s) {
     model.solveStep();
   }

An example of explicit dynamic heat propagation is presented in
``examples/heat_transfer/explicit_heat_transfer.cc``.  This example consists
of a square 2D plate of :math:`1 \text{m}^2` having an initial temperature of
:math:`100 \text{K}` everywhere but a none centered hot point maintained at
:math:`300 \text{K}`. :numref:`fig:htm:explicit:dynamic-1` presents the geometry
of this case. The material used is a linear fictitious elastic material with a
density of :math:`8940 \text{kg}/\text{m}^3`, a conductivity of
:math:`401 \text{W}/\text{m}/\text{K}` and a specific heat capacity of
:math:`385 \text{J}/\text{K}/\text{kg}`. The time step used is
:math:`0.12 \text{s}`.

.. _fig:htm:explicit:dynamic-1:
.. figure:: figures/hot-point-1.png
   :align: center

   Initial temperature field

.. _fig:htm:explicit:dynamic-2:
.. figure:: figures/hot-point-2.png
   :align: center

   Temperature field after 15000 time steps = 30 minutes. The lines represent iso-surfaces.
