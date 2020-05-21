Solid Mechanics Model
=====================

The solid mechanics model is a specific implementation of the ``Model``
interface dedicated to handle the equations of motion or equations of
equilibrium. The model is created for a given mesh. It will create its own
``FEEngine`` object to compute the interpolation, gradient, integration and
assembly operations. A :cpp:class:`SolidMechanicsModel
<akantu::SolidMechanicsModel>` object can simply be created like this::

  SolidMechanicsModel model(mesh);

where ``mesh`` is the mesh for which the equations are to be
solved. A second parameter called ``spatial_dimension`` can be
added after ``mesh`` if the spatial dimension of the problem is
different than that of the mesh.

This model contains at least the following six ``Arrays``:

blocked_dofs
    contains a Boolean value for each degree of freedom specifying whether that
    degree is blocked or not. A Dirichlet boundary condition can be prescribed
    by setting the **blocked_dofs** value of a degree of freedom to
    ``true``. A Neumann boundary condition can be applied by setting the
    **blocked_dofs** value of a degree of freedom to ``false``. The
    **displacement**, **velocity** and **acceleration** are
    computed for all degrees of freedom for which the **blocked_dofs**
    value is set to ``false``. For the remaining degrees of freedom, the imposed
    values (zero by default after initialization) are kept.

displacement
    contains the displacements of all degrees of freedom. It can be either a
    computed displacement for free degrees of freedom or an imposed displacement
    in case of blocked ones (:math:`\vec{u}` in the following).

velocity
    contains the velocities of all degrees of freedom. As **displacement**,
    it contains computed or imposed velocities depending on the nature of the
    degrees of freedom (:math:`\dot{\vec{u}}` in the following).

acceleration
    contains the accelerations of all degrees of freedom. As **displacement**,
    it contains computed or imposed accelerations depending on the nature of the
    degrees of freedom (:math:`\ddot{\vec{u}}` in the following).

external_force
    contains the external forces applied on the nodes
    (:math:`\vec{f}_{\st{ext}}` in the following).

internal_force
    contains the internal forces on the nodes (:math:`\vec{f}_{\mathrm{int}}` in
    the following).


Some examples to help to understand how to use this model will be
presented in the next sections.


Model Setup
-----------


Setting Initial Conditions
``````````````````````````

For a unique solution of the equations of motion, initial
displacements and velocities for all degrees of freedom must be
specified:

.. math::
  \vec{u}(t=0) & = \vec{u}_0\\
  \dot{\vec u}(t=0) & = \vec{v}_0

The solid mechanics model can be initialized as
follows::

  model.initFull()

This function initializes the internal arrays and sets them to
zero. Initial displacements and velocities that are not equal to zero
can be prescribed by running a loop over the total number of
nodes. Here, the initial displacement in :math:`x`-direction and the
initial velocity in :math:`y`-direction for all nodes is set to :math:`0.1` and :math:`1`,
respectively::

    auto & disp = model.getDisplacement();
    auto & velo = model.getVelocity();

    for (UInt node = 0; node < mesh.getNbNodes(); ++node) {
        disp(node, 0) = 0.1;
        velo(node, 1) = 1.;
    }

Setting Boundary Conditions
```````````````````````````

This section explains how to impose Dirichlet or Neumann boundary
conditions. A Dirichlet boundary condition specifies the values that
the displacement needs to take for every point :math:`x` at the boundary
(:math:`\Gamma_u`) of the problem domain (:numref:`fig:smm:boundaries`):

.. math::
  \vec{u} = \bar{\vec u} \quad \forall \vec{x}\in \Gamma_{u}


A Neumann boundary condition imposes the value of the gradient of the
solution at the boundary :math:`\Gamma_t` of the problem domain
(:numref:`fig:smm:boundaries`):

.. math::
    \vec{t} = \mat{\sigma} \vec{n} = \bar{\vec t} \quad
    \forall \vec{x}\in \Gamma_{t}

.. _fig:smm:boundaries:
.. figure:: figures/problem_domain.svg
            :align: center
            :width: 75%

            Problem domain :math:`\Omega` with boundary in three dimensions. The
            Dirchelet and the Neumann regions of the boundary are denoted with
            :math:`\Gamma_u` and :math:`\Gamma_t`, respecitvely.

Different ways of imposing these boundary conditions exist. A basic
way is to loop over nodes or elements at the boundary and apply local
values. A more advanced method consists of using the notion of the
boundary of the mesh. In the following both ways are presented.

Starting with the basic approach, as mentioned, the Dirichlet boundary
conditions can be applied by looping over the nodes and assigning the
required values. :numref:`fig:smm:dirichlet_bc` shows a beam with a
fixed support on the left side. On the right end of the beam, a load
is applied. At the fixed support, the displacement has a given
value. For this example, the displacements in both the :math:`x` and the
:math:`y`-direction are set to zero. Implementing this displacement boundary
condition is similar to the implementation of initial displacement
conditions described above. However, in order to impose a displacement
boundary condition for all time steps, the corresponding nodes need to
be marked as boundary nodes using the function ``blocked``. While,
in order to impose a load on the right side, the nodes are not marked.
The detail codes are shown as follows::

   auto & blocked = model.getBlockedDOFs();
   const auto & pos = mesh.getNodes();

   UInt nb_nodes = mesh.getNbNodes();

   for (UInt node = 0; node < nb_nodes; ++node) {
     if(Math::are_float_equal(pos(node, _x), 0)) {
       blocked(node, _x) = true; // block dof in x-direction
       blocked(node, _y) = true; // block dof in y-direction
       disp(node, _x) = 0.; // fixed displacement in x-direction
       disp(node, _y) = 0.; // fixed displacement in y-direction
     } else if (Math::are_float_equal(pos(node, _y), 0)) {
       blocked(node, _x) = false; // unblock dof in x-direction
       forces(node, _x) = 10.;    // force in x-direction
     }
   }


.. _fig:smm:dirichlet_bc:
.. figure:: figures/dirichlet.svg
            :align: center

            Beam with fixed support and load.


For the more advanced approach, one needs the notion of a boundary in
the mesh. Therefore, the boundary should be created before boundary
condition functors can be applied. Generally the boundary can be
specified from the mesh file or the geometry.  For the first case, the
function ``createGroupsFromMeshData`` is called.  This function
can read any types of mesh data which are provided in the mesh
file. If the mesh file is created with Gmsh, the function takes one
input strings which is either ``tag_0``, ``tag_1`` or
``physical_names``. The first two tags are assigned by Gmsh to
each element which shows the physical group that they belong to. In
Gmsh, it is also possible to consider strings for different groups of
elements. These elements can be separated by giving a string
``physical_names`` to the function
``createGroupsFromMeshData``::

    mesh.createGroupsFromMeshData<std::string>("physical_names").

Boundary conditions support can also be created from the geometry by calling
``createBoundaryGroupFromGeometry``. This function gathers all the elements on
the boundary of the geometry.

To apply the required boundary conditions, the function ``applyBC`` needs to be
called on a :cpp:class:`SolidMechanicsModel <akantu::SolidMechanicsModel>`. This
function gets a Dirichlet or Neumann functor and a string which specifies the
desired boundary on which the boundary conditions is to be applied. The functors
specify the type of conditions to apply. Three built-in functors for Dirichlet
exist: ``FlagOnly, FixedValue,`` and ``IncrementValue``. The functor
``FlagOnly`` is used if a point is fixed in a given direction. Therefore, the
input parameter to this functor is only the fixed direction. The ``FixedValue``
functor is used when a displacement value is applied in a fixed direction. The
``IncrementValue`` applies an increment to the displacement in a given
direction. The following code shows the utilization of three functors for the
top, bottom and side surface of the mesh which were already defined in the Gmsh
file::

   model.applyBC(BC::Dirichlet::FixedValue(13.0, _y), "Top");
   model.applyBC(BC::Dirichlet::FlagOnly(_x), "Bottom");
   model.applyBC(BC::Dirichlet::IncrementValue(13.0, _x), "Side");

To apply a Neumann boundary condition, the applied traction or stress should be
specified before. In case of specifying the traction on the surface, the functor
``FromTraction`` of Neumann boundary conditions is called. Otherwise, the
functor ``FromStress`` should be called which gets the stress tensor as an input
parameter::

   Vector<Real> surface_traction = {0., 0., 1.};
   auto surface_stress(3, 3) = Matrix<Real>::eye(3);

   model.applyBC(BC::Neumann::FromTraction(surface_traction), "Bottom");
   model.applyBC(BC::Neumann::FromStress(surface_stress), "Top");

If the boundary conditions need to be removed during the simulation, a
functor is called from the Neumann boundary condition to free those
boundary conditions from the desired boundary::

    model.applyBC(BC::Neumann::FreeBoundary(), "Side");

User specified functors can also be implemented.  A full example for
setting both initial and boundary conditions can be found in
``examples/boundary_conditions.cc``.  The problem solved
in this example is shown in Fig.~\ref{fig:smm:bc_and_ic}. It consists
of a plate that is fixed with movable supports on the left and bottom
side. On the right side, a traction, which increases linearly with the
number of time steps, is applied. The initial displacement and
velocity in :math:`x`-direction at all free nodes is zero and two
respectively.

.. _fig:smm:bc_and_ic:
.. figure:: figures/bc_and_ic_example.svg
            :align: center
            :width: 75%

            Plate on movable supports.

..
   \begin{figure}[!htb]
     \centering
     \includegraphics[scale=0.8]{figures/bc_and_ic_example}
     \caption{Plate on movable supports.\label{fig:smm:bc_and_ic}}
   \end{figure}

As it is mentioned in Section \ref{sect:common:groups}, node and
element groups can be used to assign the boundary conditions. A
generic example is given below with a Dirichlet boundary condition::

     // create a node group
     NodeGroup & node_group = mesh.createNodeGroup("nodes_fix");

     /* fill the node group with the nodes you want */

     // create an element group using the existing node group
     mesh.createElementGroupFromNodeGroup("el_fix",
                                          "nodes_fix",
                                          spatial_dimension-1);

     // boundary condition can be applied using the element group name
     model.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "el_fix");

Material Selector
`````````````````

If the user wants to assign different materials to different
finite elements groups in \akantu, a material selector has to be
used. By default, \akantu assigns the first valid material in the
material file to all elements present in the model (regular continuum
materials are assigned to the regular elements and cohesive materials
are assigned to cohesive elements or element facets).

To assign different materials to specific elements, mesh data
information such as tag information or specified physical names can be
used. ``MeshDataMaterialSelector`` class uses this information to
assign different materials. With the proper physical name or tag name
and index, different materials can be assigned as demonstrated in the
examples below::

     auto mat_selector =
        std::make_shared<MeshDataMaterialSelector<std::string>>("physical_names",
                                                                model);
     model.setMaterialSelector(mat_selector);

In this example the physical names specified in a GMSH geometry file will by
used to match the material names in the input file.

Another example would be to use the first (``tag_0``) or the second
(``tag_1``) tag associated to each elements in the mesh::

     auto mat_selector = std::make_shared<MeshDataMaterialSelector<UInt>>(
         "tag_1", model, first_index);
     model.setMaterialSelector(*mat_selector);

where ``first_index`` (default is 1) is the value of ``tag_1`` that will
be associated to the first material in the material input file. The following
values of the tag will be associated with the following materials.

There are four different material selectors pre-defined in
\akantu. ``MaterialSelector`` and ``DefaultMaterialSelector`` is
used to assign a material to regular elements by default. For the
regular elements, as in the example above,
``MeshDataMaterialSelector`` can be used to assign different
materials to different elements.

Apart from the \akantu's default material selectors, users can always
develop their own classes in the main code to tackle various
multi-material assignment situations.

Static Analysis
---------------

The :cpp:class:`SolidMechanicsModel <akantu::SolidMechanicsModel>` class can
handle different analysis methods, the first one being presented is the static
case. In this case, the equation to solve is

.. math::
     \mat{K} \vec{u} = \vec{f}_{\mathrm{ext}}

where :math:`\mat{K}` is the global stiffness matrix, :math:`\vec{u}` the
displacement vector and :math:`\vec{f}_{\st{ext}}` the vector of external
forces applied to the system.

To solve such a problem, the static solver of the
:cpp:class:`SolidMechanicsModel <akantu::SolidMechanicsModel>` object is used.
First, a model has to be created and initialized. To create the model, a mesh
(which can be read from a file) is needed, as explained in
Section~\ref{sect:common:mesh}. Once an instance of a
:cpp:class:`SolidMechanicsModel <akantu::SolidMechanicsModel>` is obtained, the
easiest way to initialize it is to use the ``initFull`` method by giving the
``SolidMechanicsModelOptions``. These options specify the type of analysis to be
performed and whether the materials should be initialized with ``initMaterials``
or not::

   SolidMechanicsModel model(mesh);
   model.initFull(_analysis_method = _static);

Here, a static analysis is chosen by passing the argument
``_static`` to the method. By default, the Boolean for no
initialization of the materials is set to false, so that they are
initialized during the ``initFull``. The method ``initFull``
also initializes all appropriate vectors to zero.  Once the model is
created and initialized, the boundary conditions can be set as
explained in Section~\ref{sect:smm:boundary}.  Boundary conditions
will prescribe the external forces for some free degrees of freedom
:math:`\vec{f}_{\st{ext}}` and displacements for some others.  At this point
of the analysis, the function
``solveStep``\index{SolidMechanicsModel!solveStep} can be called::

   model.solveStep<_scm_newton_raphson_tangent_modified,
                   SolveConvergenceCriteria::_residual>(1e-4, 1);

This function is templated by the solving method and the convergence criterion
and takes two arguments: the tolerance and the maximum number of iterations (100
by default), which are :math:`10^{-4}` and :math:`1` for this example. The
modified Newton-Raphson method is chosen to solve the system. In this method,
the equilibrium equation (\ref{eqn:smm:static}) is modified in order to apply a
Newton-Raphson convergence algorithm:

.. math::
     \mat{K}^{i+1}\delta\vec{u}^{i+1} &= \vec{r} \\
     &= \vec{f}_{\st{ext}} -\vec{f}_{\st{int}}\\
     &= \vec{f}_{\st{ext}} - \mat{K}^{i} \vec{u}^{i}\\
     \vec{u}^{i+1} &= \vec{u}^{i} + \delta\vec{u}^{i+1}~,

where :math:`\delta\vec{u}` is the increment of displacement to be added from
one iteration to the other, and :math:`i` is the Newton-Raphson iteration
counter. By invoking the ``solveStep`` method in the first step, the global
stiffness matrix :math:`\mat{K}` from Equation~(\ref{eqn:smm:static}) is
automatically assembled. A Newton-Raphson iteration is subsequently started,
:math:`\mat{K}` is updated according to the displacement computed at the
previous iteration and one loops until the forces are balanced
(``SolveConvergenceCriteria::_residual``), i.e. :math:`||\vec{r}|| <
\texttt{SolveConvergenceCriteria::_residual}`. One can also iterate until the
increment of displacement is zero (``SolveConvergenceCriteria::_increment``)
which also means that the equilibrium is found. For a linear elastic problem,
the solution is obtained in one iteration and therefore the maximum number of
iterations can be set to one. But for a non-linear case, one needs to iterate as
long as the norm of the residual exceeds the tolerance threshold and therefore
the maximum number of iterations has to be higher, e.g. :math:`100`::

    model.solveStep<_scm_newton_raphson_tangent_modified,
                    SolveConvergenceCriteria::_residual>(1e-4, 100)

At the end of the analysis, the final solution is stored in the
**displacement** vector.  A full example of how to solve a static
problem is presented in the code ``examples/static/static.cc``.
This example is composed of a 2D plate of steel, blocked with rollers
on the left and bottom sides as shown in :numref:`fig:smm:static`.
The nodes from the right side of the sample are displaced by :math:`0.01\%`
of the length of the plate.

.. _fig:smm:static:
.. figure:: figures/static.svg
            :align: center
            :width: 75%

            Numerical setup.

The results of this analysis is depicted in
:numref:`fig:smm:implicit:static_solution`.

..
   \begin{figure}[!htb]
     \centering
     \includegraphics[width=.7\linewidth]{figures/static_analysis}
     \caption{Solution of the static analysis. Left: the initial
   condition, right: the solution (deformation magnified 50 times)}
     \label{fig:smm:implicit:static_solution}
   \end{figure}

.. _fig:smm:implicit:static_solution:
.. figure:: figures/static_analysis.png
            :align: center
            :width: 75%

            Solution of the static analysis. Left: the initial condition, right:
            the solution (deformation magnified 50 times).

Static implicit analysis with dynamic insertion of cohesive elements
````````````````````````````````````````````````````````````````````


In order to solve problems with the extrinsic cohesive method in the
static implicit solution scheme, the function ``solveStepCohesive``
has to be used::

   model.solveStepCohesive<_scm_newton_raphson_tangent,
                           SolveConvergenceCriteria::_increment>(1e-13, error,
                                                                 25, false, 1e5,
                                                                 true);


in which the arguments are: tolerance, error, max_iteration,
load_reduction, tol_increase_factor, do_not_factorize.  This
function, first applies the Newton-Raphson procedure to solve the
problem.  Then, it calls the method ``checkCohesiveStress`` to
check if cohesive elements have to be inserted.  Since the approach is
implicit, only one element is added, the most stressed one (see
Section \ref{extrinsic_insertion}).  After insertion, the
Newton-Raphson procedure is applied again to solve the same
incremental loading step, with the new inserted cohesive element.  The
procedure loops in this way since no new cohesive elements have to be
inserted.  At that point, the solution is saved, and the simulation
can advance to the next incremental loading step.  In case the
convergence is not reached, the obtained solution is not saved and the
simulation return to the main file with the error given by the
solution saved in the argument of the function *error*.  In this
way, the user can intervene in the simulation in order to find anyhow
convergence.  A possibility is, for instance, to reduce the last
incremental loading step.  The variable *load_reduction* can be
used to identify if the load has been already reduced or not.  At the
same time, with the variable *tol_increase_factor* it is
possible to increase the tolerance by a factor defined by the user in
the main file, in order to accept a solution even with an error bigger
than the tolerance set at the beginning.  It is possible to increase
the tolerance only in the phase of loading reduction, i.e., when
load_reduction = true.  A not converged solution is never saved.  In
case the convergence is not reached even after the loading reduction
procedure, the displacement field is not updated and remains the one
of the last converged incremental steps.  Also, cohesive elements are
inserted only if convergence is reached.  An example of the extrinsic
cohesive method in the static implicit solution scheme is presented in
``examples/cohesive_element/cohesive_extrinsic_implicit``.

Dynamic Methods
---------------

Different ways to solve the equations of motion are implemented in the
solid mechanics model.  The complete equations that should be solved
are:

.. _eqn:equation-motion:
.. math::
   \mat{M}\ddot{\vec{u}} + \mat{C}\dot{\vec{u}} + \mat{K}\vec{u} =
   \vec{f}_{\mathrm{ext}}

where :math:`\mat{M}`, :math:`\mat{C}` and :math:`\mat{K}` are the mass,
damping and stiffness matrices, respectively.

In the previous section, it has already been discussed how to solve this
equation in the static case, where :math:`\ddot{\vec{u}} = \dot{\vec{u}} = 0`.  Here
the method to solve this equation in the general case will be presented.  For
this purpose, a time discretization has to be specified.  The most common
discretization method in solid mechanics is the Newmark-:math:`\beta` method, which is
also the default in Akantu.

For the Newmark-:math:`\beta` method, (:numref:`eqn:equation-motion`) becomes a
system of three equations (see \cite{curnier92a} \cite{hughes-83a} for
more details):

.. _eqn:finite-difference:
.. math::
   \mat{M} \ddot{\vec{u}}_{n+1} + \mat{C}\dot{\vec{u}}_{n+1} + \mat{K} \vec{u}_{n+1} &={\vec{f}_{\st{ext}}}_{\, n+1}
   \label{eqn:equation-motion-discret} \\
   \vec{u}_{n+1} &=\vec{u}_{n} + \left(1 - \alpha\right) \Delta t \dot{\vec{u}}_{n} +
   \alpha \Delta t \dot{\vec{u}}_{n+1} + \left(\frac{1}{2} -
   \alpha\right) \Delta t^2
   \ddot{\vec{u}}_{n} \label{eqn:finite-difference-1}\\
   \dot{\vec{u}}_{n+1} &= \dot{\vec{u}}_{n} + \left(1 - \beta\right)
   \Delta t \ddot{\vec{u}}_{n} + \beta \Delta t
   \ddot{\vec{u}}_{n+1}


In these new equations, :math:`\ddot{\vec{u}}_{n}`, :math:`\dot{\vec{u}}_{n}`
and :math:`\vec{u}_{n}` are the approximations of :math:`\ddot{\vec{u}}(t_n)`,
:math:`\dot{\vec{u}}(t_n)` and :math:`\vec{u}(t_n)`.
Equation~(\ref{eqn:equation-motion-discret}) is the equation of motion
discretized in space (finite-element discretization), and the equations above
are discretized in both space and time (Newmark discretization). The
:math:`\alpha` and :math:`\beta` parameters determine the stability and the
accuracy of the algorithm. Classical values for :math:`\alpha` and :math:`\beta`
are usually :math:`\beta = 1/2` for no numerical damping and :math:`0 < \alpha <
1/2`.

+-------------------+-------------+----------+
|:math:`\alpha`     |Method       |Type      |
|                   |(:math:`\beta|          |
|                   |= 1/2`)      |          |
+-------------------+-------------+----------+
|:math:`0`          |central      |explicit  |
|                   |difference   |          |
+-------------------+-------------+----------+
|:math:`\frac{1}{6}`|Fox-Goodwin  |implicit  |
|                   |(royal road) |          |
+-------------------+-------------+----------+
|:math:`\frac{1}{3}`|Linear       |implicit  |
|                   |acceleration |          |
+-------------------+-------------+----------+
|:math:`\frac{1}{2}`|Average      |implicit  |
|                   |acceleration |          |
|                   |(trapeziodal |          |
|                   |rule)        |          |
+-------------------+-------------+----------+


The solution of this system of equations,
(\ref{eqn:equation-motion-discret})-(\ref{eqn:finite-difference-2}) is
split into a predictor and a corrector system of equations.  Moreover,
in the case of a non-linear equations, an iterative algorithm such as
the Newton-Raphson method is applied. The system of equations can be
written as:

- *Predictor:*
.. math::
    \vec{u}_{n+1}^{0} &= \vec{u}_{n} + \Delta t
    \dot{\vec{u}}_{n} + \frac{\Delta t^2}{2} \ddot{\vec{u}}_{n} \\
    \dot{\vec{u}}_{n+1}^{0} &= \dot{\vec{u}}_{n} + \Delta t
    \ddot{\vec{u}}_{n} \\
    \ddot{\vec{u}}_{n+1}^{0} &= \ddot{\vec{u}}_{n}
- *Solve:*
.. math::

     \left(c \mat{M} + d \mat{C} + e \mat{K}_{n+1}^i\right)
     \vec{w} = {\vec{f}_{\st{ext}}}_{\,n+1} - {\vec{f}_{\st{int}}}_{\,n+1}^i -
     \mat{C} \dot{\vec{u}}_{n+1}^i - \mat{M} \ddot{\vec{u}}_{n+1}^i = \vec{r}_{n+1}^i
- *Corrector:*
.. math::
     \ddot{\vec{u}}_{n+1}^{i+1} &= \ddot{\vec{u}}_{n+1}^{i} +c \vec{w} \\
     \dot{\vec{u}}_{n+1}^{i+1} &= \dot{\vec{u}}_{n+1}^{i} + d\vec{w} \\
     \vec{u}_{n+1}^{i+1} &= \vec{u}_{n+1}^{i} + e \vec{w}

where :math:`i` is the Newton-Raphson iteration counter and :math:`c`, :math:`d` and :math:`e`
are parameters depending on the method used to solve the equations

+--------------------+----------------------------+------------------------+----------------------+----------------------+
|                    |:math:`\vec{w}`             |:math:`e`               |:math:`d`             |:math:`c`             |
|                    |                            |                        |                      |                      |
|                    |                            |                        |                      |                      |
+--------------------+----------------------------+------------------------+----------------------+----------------------+
|in acceleration     |:math:`\delta\ddot{\vec{u}}`|:math:`\alpha\beta\Delta|:math:`\beta\Delta    |:math:`1`             |
|                    |                            |t^2`                    |t`                    |                      |
|                    |                            |                        |                      |                      |
+--------------------+----------------------------+------------------------+----------------------+----------------------+
|in velocity         |:math:`\delta\dot{\vec{u}}` |:math:`\alpha\Delta t`  |:math:`1`             |:math:`\frac{1}{\beta |
|                    |                            |                        |                      |\Delta t}`            |
|                    |                            |                        |                      |                      |
+--------------------+----------------------------+------------------------+----------------------+----------------------+
|in displacement     |:math:`\delta\vec{u}`       |:math:`1`               |:math:`\frac{1}{\alpha|:math:`\frac{1}{\alpha|
|                    |                            |                        |\Delta t}`            |\beta \Delta t^2}`    |
|                    |                            |                        |                      |                      |
+--------------------+----------------------------+------------------------+----------------------+----------------------+


..
   \begin{center}
     \begin{tabular}{lcccc}
       \toprule
       & :math:`\vec{w}` & :math:`e` & :math:`d` & :math:`c`\\
       \midrule
       in acceleration &:math:` \delta\ddot{\vec{u}}` & :math:`\alpha \beta\Delta t^2` &:math:`\beta \Delta t` &:math:`1`\\
       in velocity & :math:` \delta\dot{\vec{u}}`& :math:`\alpha\Delta t` & :math:`1` & :math:`\frac{1}{\beta \Delta t}`\\
       in displacement &:math:`\delta\vec{u}` & :math:` 1` & :math:`\frac{1}{\alpha \Delta t}` & :math:`\frac{1}{\alpha \beta \Delta t^2}`\\
       \bottomrule
     \end{tabular}
   \end{center}

..
   % \note{If you want to use the implicit solver \akantu should be compiled at
   % least with one sparse matrix solver such as Mumps\cite{mumps}.}


Implicit Time Integration
`````````````````````````

To solve a problem with an implicit time integration scheme, first a
:cpp:class:`SolidMechanicsModel <akantu::SolidMechanicsModel>` object has to be
created and initialized. Then the initial and boundary conditions have to be
set. Everything is similar to the example in the static case
(Section~\ref{sect:smm:static}), however, in this case the implicit dynamic
scheme is selected at the initialization of the model::

   SolidMechanicsModel model(mesh);
   model.initFull(_analysis_method = _implicit_dynamic);

Because a dynamic simulation is conducted, an integration time step
:math:`\Delta t` has to be specified. In the case of implicit simulations,
\akantu implements a trapezoidal rule by default.  That is to say
:math:`\alpha = 1/2` and :math:`\beta = 1/2` which is unconditionally
stable. Therefore the value of the time step can be chosen arbitrarily
within reason::

   model.setTimeStep(time_step);

Since the system has to be solved for a given amount of time steps, the
method ``solveStep()``, (which has already been used in the static
example in Section~\ref{sect:smm:static}), is called inside a time
loop::

   /// time loop
   Real time = 0.;

   auto & solver = model.getNonLinearSolver();
   solver.set("max_iterations", 100);
   solver.set("threshold", 1e-12);
   solver.set("convergence_type", SolveConvergenceCriteria::_solution);

   for (UInt s = 1; time <max_time; ++s, time += time_step) {
       model.solveStep();
   }

An example of solid mechanics with an implicit time integration scheme is
presented in ``examples/implicit/implicit_dynamic.cc``. This example consists of
a 3D beam of
:math:`10\mathrm{m}\times1\mathrm{m}\times1\mathrm{m}` blocked
on one side and is on a roller on the other side. A constant force of
:math:`5\mathrm{kN}` is applied in its middle.
:numref:`fig:smm:implicit:dynamic` presents the geometry of this case. The
material used is a fictitious linear elastic material with a density of
:math:`1000 \mathrm{kg/m}^3`, a Young's Modulus of
:math:`\SI{120}{\mega\pascal}` and Poisson's ratio of :math:`0.3`. These values
were chosen to simplify the analytical solution.

An approximation of the dynamic response of the middle point of the
beam is given by:

.. math::

    u\left(\frac{L}{2}, t\right)
    = \frac{1}{\pi^4} \left(1 - cos\left(\pi^2 t\right) +
    \frac{1}{81}\left(1 - cos\left(3^2 \pi^2 t\right)\right) +
    \frac{1}{625}\left(1 - cos\left(5^2 \pi^2 t\right)\right)\right)

.. _fig:smm:implicit:dynamic:
.. figure:: figures/implicit_dynamic.svg
            :align: center
            :width: 75%

            Numerical setup.

..
   \begin{figure}[!htb]
     \centering
     \includegraphics[scale=.6]{figures/implicit_dynamic}
     \caption{Numerical setup}
     \label{fig:smm:implicit:dynamic}
   \end{figure}

Figure :numref:`fig:smm:implicit:dynamic_solution` presents the deformed
beam at 3 different times during the simulation: time steps 0, 1000 and
2000.

.. _fig:smm:implicity:dynamic_solution:
.. figure:: figures/dynamic_analysis.png
            :align: center
            :width: 60%

            Deformed beam at three different times (displacement :math:`\times
            10`).
..
   \begin{figure}[!htb]
     \centering
     \setlength{\unitlength}{0.1\textwidth}
     \begin{tikzpicture}
       \node[above right] (img) at (0,0)
       {\includegraphics[width=.6\linewidth]{figures/dynamic_analysis}};
       \node[left] at (0pt,20pt) {:math:`0`}; \node[left] at (0pt,60pt) {:math:`1000`};
       \node[left] at (0pt,100pt) {:math:`2000`};
     \end{tikzpicture}

     \caption{Deformed beam at 3 different times (displacement are
       magnified by a factor 10).}
     \label{fig:smm:implicit:dynamic_solution}
   \end{figure}

Explicit Time Integration
`````````````````````````

The explicit dynamic time integration scheme is based on the
Newmark-:math:`\beta` scheme with :math:`\alpha=0` (see equations
\ref{eqn:equation-motion-discret}-\ref{eqn:finite-difference-2}).  In
\akantu, :math:`\beta` is defaults to :math:`\beta=1/2`, see section
\ref{sect:smm:Dynamic_methods}.

The initialization of the simulation is similar to the static and
implicit dynamic version.  The model is created from the
:cpp:class:`SolidMechanicsModel <akantu::SolidMechanicsModel>` class.  In the initialization, the explicit
scheme is selected using the ``_explicit_lumped_mass`` constant::

   SolidMechanicsModel model(mesh);
   model.initFull(_analysis_method = _explicit_lumped_mass);


.. note::
    Writing ``model.initFull()`` or ``model.initFull();`` is
    equivalent to use the ``_explicit_lumped_mass`` keyword, as this
    is the default case.

The explicit time integration scheme implemented in \akantu uses a
lumped mass matrix :math:`\mat{M}` (reducing the computational cost). This
matrix is assembled by distributing the mass of each element onto its
nodes. The resulting :math:`\mat{M}` is therefore a diagonal matrix stored
in the **mass** vector of the model.

The explicit integration scheme is conditionally stable. The time step
has to be smaller than the stable time step which is obtained in
Akantu as follows::

   critical_time_step = model.getStableTimeStep();

The stable time  step corresponds to the time the fastest wave (the compressive
wave) needs to travel the characteristic length of the mesh:

.. math::
   \Delta t_{\st{crit}} = \frac{\Delta x}{c}

where :math:`\Delta x` is a characteristic length (\eg the inradius in the case
of linear triangle element) and :math:`c` is the celerity of the fastest wave in
the material. It is generally the compressive wave of celerity :math:`c =
\sqrt{\frac{2 \mu + \lambda}{\rho}}`, :math:`\mu` and :math:`\lambda` are the
first and second Lame's coefficients and :math:`\rho` is the density. However,
it is recommended to impose a time step that is smaller than the stable time
step, for instance, by multiplying the stable time step by a safety factor
smaller than one::

   const Real safety_time_factor = 0.8;
   Real applied_time_step = critical_time_step * safety_time_factor;
   model.setTimeStep(applied_time_step);

The initial displacement and velocity fields are, by default, equal to zero if
not given specifically by the user (see \ref{sect:smm:initial_condition}).

Like in implicit dynamics, a time loop is used in which the
displacement, velocity and acceleration fields are updated at each
time step. The values of these fields are obtained from the
Newmark:math:`-\beta` equations with :math:`\beta=1/2` and :math:`\alpha=0`. In \akantu
these computations at each time step are invoked by calling the
function ``solveStep``::

   for (UInt s = 1; (s-1)*applied_time_step < total_time; ++s) {
     model.solveStep();
   }

The method ``solveStep`` wraps the four following functions:

- ``model.explicitPred()`` allows to compute the displacement
     field at :math:`t+1` and a part of the velocity field at :math:`t+1`, denoted by
     :math:`\vec{\dot{u}^{\st{p}}}_{n+1}`, which will be used later in the method
     ``model.explicitCorr()``. The equations are:

     .. math::
        \vec{u}_{n+1} &= \vec{u}_{n} + \Delta t
        \vec{\dot{u}}_{n} + \frac{\Delta t^2}{2} \vec{\ddot{u}}_{n}\\
        \vec{\dot{u}^{\st{p}}}_{n+1} &= \vec{\dot{u}}_{n} + \Delta t
        \vec{\ddot{u}}_{n}

- ``model.updateResidual()`` and ``model.updateAcceleration()`` compute the acceleration increment
     :math:`\delta \vec{\ddot{u}}`:

     .. math::
        \left(\mat{M} + \frac{1}{2} \Delta t \mat{C}\right)
        \delta \vec{\ddot{u}} = \vec{f_{\st{ext}}} - \vec{f}_{\st{int}\, n+1}
        - \mat{C} \vec{\dot{u}^{\st{p}}}_{n+1} - \mat{M} \vec{\ddot{u}}_{n}

     The internal force :math:`\vec{f}_{\st{int}\, n+1}` is computed from the
       displacement :math:`\vec{u}_{n+1}` based on the constitutive law.

- ``model.explicitCorr()`` computes the velocity and
     acceleration fields at :math:`t+1`:

     .. math::
        \vec{\dot{u}}_{n+1} &= \vec{\dot{u}^{\st{p}}}_{n+1} + \frac{\Delta t}{2}
        \delta \vec{\ddot{u}} \\ \vec{\ddot{u}}_{n+1} &=
        \vec{\ddot{u}}_{n} + \delta \vec{\ddot{u}}

The use of an explicit time integration scheme is illustrated by the example:
``examples/explicit/explicit_dynamic.cc``. This example models the propagation
of a wave in a steel beam. The beam and the applied displacement in the
:math:`x` direction are shown in :numref:`fig:smm:explicit`.


.. _fig:smm:explicit:
.. figure:: figures/explicit.svg
            :align: center
            :width: 90%

            Numerical setup.

..
   \begin{figure}[!htb] \centering
     \begin{tikzpicture}
       \coordinate (c) at (0,2);
       \draw[shift={(c)},thick, color=blue] plot [id=x, domain=-5:5, samples=50] ({\x, {(40 * sin(0.1*pi*3*\x) * exp(- (0.1*pi*3*\x)*(0.1*pi*3*\x) / 4))}});
       \draw[shift={(c)},-latex] (-6,0) -- (6,0) node[right, below] {:math:`x`};
       \draw[shift={(c)},-latex] (0,-0.7) -- (0,1) node[right] {:math:`u`};
       \draw[shift={(c)}] (-0.1,0.6) node[left] {:math:`A`}-- (1.5,0.6);

       \coordinate (l) at (0,0.6);
       \draw[shift={(0,-0.7)}] (-5, 0) -- (5,0) -- (5, 1) -- (-5, 1) -- cycle;
       \draw[shift={(l)}, latex-latex] (-5,0)-- (5,0) node [midway, above] {:math:`L`};
       \draw[shift={(l)}] (5,0.2)-- (5,-0.2);
       \draw[shift={(l)}] (-5,0.2)-- (-5,-0.2);

       \coordinate (h) at (5.3,-0.7);
       \draw[shift={(h)}, latex-latex] (0,0)-- (0,1) node [midway, right] {:math:`h`};
       \draw[shift={(h)}] (-0.2,1)-- (0.2,1);
       \draw[shift={(h)}] (-0.2,0)-- (0.2,0);
     \end{tikzpicture}

     \caption{Numerical setup \label{fig:smm:explicit}}
   \end{figure}

The length and height of the beam are :math:`L={10}\textrm{m}` and :math:`h =
{1}\textrm{m}`, respectively. The material is linear elastic, homogeneous and
isotropic (density: \SI{7800}{\kilo\gram\per\cubic\metre}, Young's modulus:
\SI{210}{\giga\pascal} and Poisson's ratio: :math:`0.3`). The imposed
displacement follow a Gaussian function with a maximum amplitude of :math:`A =
{0.01}\textrm{m}`. The potential, kinetic and total energies are computed. The
safety factor is equal to :math:`0.8`.

   \input{manual-constitutive-laws}

Adding a New Constitutive Law
-----------------------------

There are several constitutive laws in \akantu as described in the
previous Section~\ref{sect:smm:CL}. It is also possible to use a
user-defined material for the simulation. These materials are referred
to as local materials since they are local to the example of the user
and not part of the \akantu library.  To define a new local material,
two files (``material_XXX.hh`` and ``material_XXX.cc``) have
to be provided where ``XXX`` is the name of the new material. The
header file ``material_XXX.hh`` defines the interface of your
custom material. Its implementation is provided in the
``material_XXX.cc``. The new law must inherit from the
``Material`` class or any other existing material class. It is
therefore necessary to include the interface of the parent material
in the header file of your local material and indicate the inheritance
in the declaration of the class::

   /* ---------------------------------------------------------------------- */
   #include "material.hh"
   /* ---------------------------------------------------------------------- */

   #ifndef __AKANTU_MATERIAL_XXX_HH__
   #define __AKANTU_MATERIAL_XXX_HH__

   namespace akantu {

   class MaterialXXX : public Material {

   /// declare here the interface of your material

   };

In the header file the user also needs to declare all the members of the new
material. These include the parameters that a read from the
material input file, as well as any other material parameters that will be
computed during the simulation and internal variables.


In the following the example of adding a new damage material will be
presented. In this case the parameters in the material will consist of the
Young's modulus, the Poisson coefficient, the resistance to damage and the
damage threshold. The material will then from these values compute its Lamé
coefficients and its bulk modulus. Furthermore, the user has to add a new
internal variable ``damage`` in order to store the amount of damage at each
quadrature point in each step of the simulation. For this specific material the
member declaration inside the class will look as follows::

   class LocalMaterialDamage : public Material {

   /// declare constructors/destructors here

   /// declare methods and accessors here

     /* -------------------------------------------------------------------- */
     /* Class Members                                                        */
     /* -------------------------------------------------------------------- */

     AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Damage, damage, Real);
   private:

     /// the young modulus
     Real E;

     /// Poisson coefficient
     Real nu;

     /// First Lame coefficient
     Real lambda;

     /// Second Lame coefficient (shear modulus)
     Real mu;

     /// resistance to damage
     Real Yd;

     /// damage threshold
     Real Sd;

     /// Bulk modulus
     Real kpa;

     /// damage internal variable
     InternalField<Real> damage;

   };

In order to enable to print the material parameters at any point in
the user's example file using the standard output stream by typing::

   for (UInt m = 0; m  < model.getNbMaterials(); ++m)
     std::cout << model.getMaterial(m) << std::endl;

the standard output stream operator has to be redefined. This should be done at the end of the header file::

   class LocalMaterialDamage : public Material {

     /// declare here the interace of your material

   }:
   /* ---------------------------------------------------------------------- */
   /* inline functions                                                       */
   /* ---------------------------------------------------------------------- */
   /// standard output stream operator
   inline std::ostream & operator <<(std::ostream & stream, const LocalMaterialDamage & _this)
   {
     _this.printself(stream);
     return stream;
   }

However, the user still needs to register the material parameters that
should be printed out. The registration is done during the call of the
constructor. Like all definitions the implementation of the
constructor has to be written in the ``material_XXX.cc``
file. However, the declaration has to be provided in the
``material_XXX.hh`` file::

   class LocalMaterialDamage : public Material {
     /* -------------------------------------------------------------------- */
     /* Constructors/Destructors                                             */
     /* -------------------------------------------------------------------- */
   public:

     LocalMaterialDamage(SolidMechanicsModel & model, const ID & id = "");
   };

The user can now define the implementation of the constructor in the
``material_XXX.cc`` file::

   /* ---------------------------------------------------------------------- */
   #include "local_material_damage.hh"
   #include "solid_mechanics_model.hh"

   namespace akantu {

   /* ---------------------------------------------------------------------- */
   LocalMaterialDamage::LocalMaterialDamage(SolidMechanicsModel & model,
              const ID & id)  :
     Material(model, id),
     damage("damage", *this) {
     AKANTU_DEBUG_IN();

     this->registerParam("E", E, 0., _pat_parsable, "Young's modulus");
     this->registerParam("nu", nu, 0.5, _pat_parsable, "Poisson's ratio");
     this->registerParam("lambda", lambda, _pat_readable, "First Lame coefficient");
     this->registerParam("mu", mu, _pat_readable, "Second Lame coefficient");
     this->registerParam("kapa", kpa, _pat_readable, "Bulk coefficient");
     this->registerParam("Yd", Yd,   50., _pat_parsmod);
     this->registerParam("Sd", Sd, 5000., _pat_parsmod);

     damage.initialize(1);

     AKANTU_DEBUG_OUT();
   }

During the intializer list the reference to the model and the material id are
assigned and the constructor of the internal field is called. Inside the scope
of the constructor the internal values have to be initialized and the
parameters, that should be printed out, are registered with the function:
``registerParam``::

   void registerParam(name of the parameter (key in the material file),
          member variable,
          default value (optional parameter),
          access permissions,
          description);

The available access permissions are as follows:
- ``_pat_internal``: Parameter can only be output when the material is printed.
- ``_pat_writable``: User can write into the parameter. The parameter is output when the material is printed.
- ``_pat_readable``: User can read the parameter. The parameter is output when the material is printed.
- ``_pat_modifiable``: Parameter is writable and readable.
- ``_pat_parsable``: Parameter can be parsed, *i.e.* read from the input file.
- ``_pat_parsmod``: Parameter is modifiable and parsable.

In order to implement the new constitutive law the user needs to
specify how the additional material parameters, that are not
defined in the input material file, should be calculated. Furthermore,
it has to be defined how stresses and the stable time step should be
computed for the new local material. In the case of implicit
simulations, in addition, the computation of the tangent stiffness needs
to be defined. Therefore, the user needs to redefine the following
functions of the parent material::

   void initMaterial();

   // for explicit and implicit simulations void
   computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

   // for implicit simulations
   void computeTangentStiffness(const ElementType & el_type,
              Array<Real> & tangent_matrix,
              GhostType ghost_type = _not_ghost);

   // for explicit and implicit simulations
   Real getStableTimeStep(Real h, const Element & element);

In the following a detailed description of these functions is provided:
- ``initMaterial``: This method is called after the material file is fully read
     and the elements corresponding to each material are assigned. Some of the
     frequently used constant parameters are calculated in this method. For
     example, the Lam\'{e} constants of elastic materials can be considered as
     such parameters.

- ``computeStress``: In this method, the stresses are computed based on the
     constitutive law as a function of the strains of the quadrature points. For
     example, the stresses for the elastic material are calculated based on the
     following formula:
     .. math::
        \mat{\sigma }  =\lambda\mathrm{tr}(\mat{\varepsilon})\mat{I}+2 \mu \mat{\varepsilon}

     Therefore, this method contains a loop on all quadrature points assigned to
     the material using the two macros:
     ``MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN`` and
     ``MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END``

     .. code::

       MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(element_type);

       // sigma <- f(grad_u)

       MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

     The strain vector in Akantu contains the values of :math:`\nabla \vec{u}`,
     i.e. it is really the *displacement gradient*,

- ``computeTangentStiffness``: This method is called when the tangent to the
     stress-strain curve is desired (see Fig \ref {fig:smm:AL:K}). For example,
     it is called in the implicit solver when the stiffness matrix for the
     regular elements is assembled based on the following formula:

     .. math::
        \label{eqn:smm:constitutive_elasc} \mat{K }
        =\int{\mat{B^T}\mat{D(\varepsilon)}\mat{B}}

     Therefore, in this method, the ``tangent`` matrix (\mat{D}) is
     computed for a given strain.

     The ``tangent`` matrix is a :math:`4^{th}` order tensor which is stored as
     a matrix in Voigt notation.

     .. _fig:smm:AL:K:
     .. figure:: figures/tangent.svg
                 :align: center
                 :width: 60%

                 Tangent to the stress-strain curve.

..
     \begin{figure}[!htb]
       \begin{center}
         \includegraphics[width=0.4\textwidth,keepaspectratio=true]{figures/tangent.pdf}
         \caption{Tangent to the stress-strain curve.}
         \label{fig:smm:AL:K}
       \end{center}
     \end{figure}

- ``getCelerity``: The stability criterion of the explicit integration scheme
  depend on the fastest wave celerity~\eqref{eqn:smm:explicit:stabletime}. This
  celerity depend on the material, and therefore the value of this velocity
  should be defined in this method for each new material. By default, the
  fastest wave speed is the compressive wave whose celerity can be defined in ``getPushWaveSpeed``.

Once the declaration and implementation of the new material has been
completed, this material can be used in the user's example by including the header file::

   #include "material_XXX.hh"

For existing materials, as mentioned in Section~\ref{sect:smm:CL}, by
default, the materials are initialized inside the method
``initFull``. If a local material should be used instead, the
initialization of the material has to be postponed until the local
material is registered in the model. Therefore, the model is
initialized with the boolean for skipping the material initialization
equal to true::

   /// model initialization
   model.initFull(_analysis_method = _explicit_lumped_mass);

Once the model has been initialized, the local material needs
to be registered in the model::

   model.registerNewCustomMaterials<XXX>("name_of_local_material");

Only at this point the material can be initialized::

   model.initMaterials();

A full example for adding a new damage law can be found in
``examples/new_material``.

Adding a New Non-Local Constitutive Law
```````````````````````````````````````

In order to add a new non-local material we first have to add the local
constitutive law in Akantu (see above). We can then add the non-local version
of the constitutive law by adding the two files (``material_XXX_non_local.hh``
and ``material_XXX_non_local.cc``) where ``XXX`` is the name of the
corresponding local material. The new law must inherit from the two classes,
non-local parent class, such as the ``MaterialNonLocal`` class, and from the
local version of the constitutive law, *i.e.* ``MaterialXXX``. It is therefore
necessary to include the interface of those classes in the header file of your
custom material and indicate the inheritance in the declaration of the class::

   /* ---------------------------------------------------------------------- */
   #include "material_non_local.hh" // the non-local parent
   #include "material_XXX.hh"
   /* ---------------------------------------------------------------------- */

   #ifndef __AKANTU_MATERIAL_XXX_HH__
   #define __AKANTU_MATERIAL_XXX_HH__

   namespace akantu {

   class MaterialXXXNonLocal : public MaterialXXX,
                               public MaterialNonLocal {

   /// declare here the interface of your material

   };

As members of the class we only need to add the internal fields to store the
non-local quantities, which are obtained from the averaging process::
  
   /* -------------------------------------------------------------------------- */
   /* Class members                                                              */
   /* -------------------------------------------------------------------------- */
   protected:
     InternalField<Real> grad_u_nl;

The following four functions need to be implemented in the non-local material::

     /// initialization of the material
     void initMaterial();
     /// loop over all element and invoke stress computation
     virtual void computeNonLocalStresses(GhostType ghost_type);
     /// compute stresses after local quantities have been averaged
     virtual void computeNonLocalStress(ElementType el_type, GhostType ghost_type)
     /// compute all local quantities
     void computeStress(ElementType el_type, GhostType ghost_type);

In the intialization of the non-local material we need to register the local
quantity for the averaging process. In our example the internal field
*grad_u_nl* is the non-local counterpart of the gradient of the displacement
field (*grad_u_nl*)::
  
     void MaterialXXXNonLocal::initMaterial() {
       MaterialXXX::initMaterial();
       MaterialNonLocal::initMaterial();
       /// register the non-local variable in the manager
       this->model->getNonLocalManager().registerNonLocalVariable(
         this->grad_u.getName(),
         this->grad_u_nl.getName(),
         spatial_dimension * spatial_dimension);
     }

The function to register the non-local variable takes as parameters the name of
the local internal field, the name of the non-local counterpart and the number
of components of the field we want to average. In the *computeStress* we now
need to compute all the quantities we want to average. We can then write a loop
for the stress computation in the function *computeNonLocalStresses* and then
provide the constitutive law on each integration point in the function
*computeNonLocalStress*.
