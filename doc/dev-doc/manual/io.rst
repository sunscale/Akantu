.. _sect-io:

Input/Output
============

Input file
----------

The text input file of a simulation should be precised using the method
:cpp:func:`initialize <akantu::initialize>` which will instantiate the static
:cpp:class:`Parser <akantu::Parser>` object of ``Akantu``. This section
explains how to manipulate at :cpp:class:`Parser <akantu::Parser>` objects to
input data in ``Akantu``.

Akantu Parser
~~~~~~~~~~~~~

``Akantu`` file parser has a tree organization.

- :cpp:class:`Parser <akantu::Parser>`, the root of the tree, can be accessed
   using::

     auto & parser = getStaticParser();

- :cpp:class:`ParserSection <akantu::ParserSection>`, branch of the tree,
   contains map a of sub-sections (:cpp:enum:`SectionType
   <akantu::SectionType>`, :cpp:class:`ParserSection <akantu::ParserSection>`)
   and a :cpp:class:`ParserSection * <akantu::ParserSection>` pointing to the
   parent section. The user section of the input file can directly be accessed
   by::

     const auto & usersect = getUserParser();

- :cpp:class:`ParserParameter <akantu::ParserParameter>`, the leaf of the tree,
   carries data of the input file which can be cast to the correct type with the
   assignment operator::

     Real mass = usersect.getParameter("mass");

   or used directly within an expression

Grammar
~~~~~~~

The structure of text input files consists of different sections
containing a list of parameters. As example, the file parsed in the
previous section will look like::

  user parameters [
    mass = 10.5
  ]

Basically every standard arithmetic operations can be used inside of input files
as well as the constant ``pi`` and ``e`` and the exponent operator ``^``.
Operations between :cpp:class:`ParserParameter <akantu::ParserParameter>` are
also possible with the convention that only parameters of the current and the
parent sections are available. :cpp:class:`Vector <akantu::Vector>` and
:cpp:class:`Matrix <akantu::Matrix>` can also be read according to the ``NumPy``
:cite:`numpy` writing convention (a.e. cauchy_stress_tensor = [[:math:`\sigma_{xx}`,
:math:`\sigma_{xy}`],[:math:`\sigma_{yx}`,\ :math:`\sigma_{yy}`]]). An example
illustrating how to parse the following input file can be found in
``example\io\parser\example_parser.cc``::

  user parameters [
      spatial_dimension = 2
      mesh_file     = swiss_cheese.msh
      inner_holes   = holes
      outter_crust  = crust
      lactostatic_p = 30e3
      stress        = [[lactostatic_p, 0            ],
                       [0,             lactostatic_p]]
      max_nb_iterations = 100
      precision     = 1e-9
  ]

.. _sect-io-material:

Material section
~~~~~~~~~~~~~~~~

The input file should also be used to specify material characteristics
(constitutive behavior and material properties). The dedicated material section
is then read by :cpp:func:`initFull <akantu::SolidMechanicsModel::initFull>`
method of :cpp:class:`SolidMechanicsModel <akantu::SolidMechanicsModel>` which
initializes the different materials specified with the following convention::

  material constitutive_law [
    name = value
    rho = value
    ...
  ]

where *constitutive_law* is the adopted constitutive law, followed by
the material properties listed one by line in the bracket (*e.g.*,
``name`` and density :math:``rho``. Some constitutive laws can also
have an *optional flavor*. More information can be found in sections
relative to material :ref:`sect-smm-cl`
or in Appendix :ref:`app-material-parameters`.

Output data
-----------

Generic data
~~~~~~~~~~~~

In this section, we address ways to get the internal data in human-readable
formats. The models in ``Akantu`` handle data associated to the mesh, but this
data can be split into several :cpp:class:`Arrays <akantu::Array>`. For example,
the data stored per element type in a :cpp:class:`ElementTypeMapArray
<akantu::ElementTypeMapArray>` is composed of as many :cpp:class:`Arrays
<akantu::Array>` as types in the mesh.

In order to get this data in a visualization software, the models contain a
object to dump ``VTK`` files. These files can be visualized in software such
as ``ParaView`` :cite:`paraview`, ``ViSit`` :cite:`visit` or ``Mayavi``
:cite:`mayavi`.

The internal dumper of the model can be configured to specify which data
fields are to be written. This is done with the  :cpp:func:`addDumpField <akantu::Model::addDumpField>` method. By default all the
files are generated in a folder called ``paraview/``::

  model.setBaseName("output"); // prefix for all generated files
  model.addDumpField("displacement"); model.addDumpField("stress"); ...
  model.dump()

The fields are dumped with the number of components of the memory. For example,
in 2D, the memory has :cpp:class:`Vectors <akantu::Vector>` of 2 components, or
the :math:`2^{nd}` order tensors with :math:`2\times2` components. This memory
can be dealt with :cpp:func:`addDumpFieldVector
<akantu::Model::addDumpFieldVector>` which always dumps :cpp:class:`Vectors
<akantu::Vector>` with 3 components or :cpp:func:`addDumpFieldTensor
<akantu::Model::addDumpFieldTensor>` which dumps :math:`2^{nd}` order tensors
with :math:`3\times3` components respectively. The routines :cpp:func:`addDumpFieldVector <akantu::Model::addDumpFieldVector>` and
:cpp:func:`addDumpFieldTensor <akantu::Model::addDumpFieldTensor>` were
introduced because of ``ParaView`` which mostly manipulate 3D data.

Those fields which are stored by quadrature point are modified to be seen in the
``VTK`` file as elemental data. To do this, the default is to average the
values of all the quadrature points.

The list of fields depends on the models (for :cpp:class:`SolidMechanicsModel
<akantu::SolidMechanicsModel>` see table :ref:`tab-io-smm-field-list`.

.. container::
   :name: tab-io-smm-field-list

   .. table:: List of dumpable fields for :cpp:class:`SolidMechanicsModel <akantu::SolidMechanicsModel>`.

      ====================== ============ =================
      key                    type         support
      ====================== ============ =================
      displacement           Vector<Real> nodes
      mass                   Vector<Real> nodes
      velocity               Vector<Real> nodes
      acceleration           Vector<Real> nodes
      force                  Vector<Real> nodes
      residual               Vector<Real> nodes
      increment              Vector<Real> nodes
      blocked_dofs           Vector<bool> nodes
      partitions             Real         elements
      material_index         variable     elements
      strain                 Matrix<Real> quadrature points
      Green strain           Matrix<Real> quadrature points
      principal strain       Vector<Real> quadrature points
      principal Green strain Vector<Real> quadrature points
      grad_u                 Matrix<Real> quadrature points
      stress                 Matrix<Real> quadrature points
      Von Mises stress       Real         quadrature points
      material_index         variable     quadrature points
      ====================== ============ =================

Cohesive elements’ data
~~~~~~~~~~~~~~~~~~~~~~~

Cohesive elements and their relative data can be easily dumped thanks to
a specific dumper contained in :cpp:class:`SolidMechanicsModelCohesive <akantu::SolidMechanicsModelCohesive>`. In order to use it, one has
just to add the string ``cohesive elements`` when calling each method already
illustrated. Here is an example on how to dump displacement and damage::

  model.addDumpFieldVectorToDumper("cohesive elements", "displacement");
  model.addDumpFieldToDumper("cohesive elements", "damage");
  model.dump("cohesive elements");


Fragmentation data
^^^^^^^^^^^^^^^^^^

Whenever the :cpp:class:`SolidMechanicsModelCohesive
<akantu::SolidMechanicsModelCohesive>` is used, it is possible to dump
additional data about the fragments that get formed in the simulation both in
serial and parallel. This task is carried out by the
:cpp:class:`FragmentManager <akantu::FragmentManager>` class, that takes care of
computing the following quantities for each fragment:

-  index;

-  mass;

-  moments of inertia;

-  velocity;

-  number of elements.

These computations can be realized at once by calling the function
:cpp:class:`computeAllData <akantu::FragmentManager::computeAllData>`, or
individually by calling the other public functions of the class. The data can be
dumped to be visualized in ``ParaView``, or can be accessed within the
simulation. An example of usage is:

At the end of this example the velocities of the fragments are accessed with a
reference to a :cpp:class:`const Array\<Real\> <akantu::Array>`. The size of this
array is the number of fragments, and its number of components is the spatial
dimension in this case.

Advanced dumping
~~~~~~~~~~~~~~~~

Arbitrary fields
^^^^^^^^^^^^^^^^

In addition to the predetermined fields from the models and materials,
the user can add any data to a dumper as long as the support is the
same. That is to say data that have the size of the full mesh on if the
dumper is dumping the mesh, or of the size of an element group if it is
a filtered dumper.

For this the easiest is to use the “external” fields register functions

The simple case force nodal and elemental data are to pass directly the
data container itself if it as the good size.

-  For nodal fields:

   It is assumed that the array as the same size as the number of nodes
   in the mesh

-  For elemental fields:

   It is assumed that the arrays in the map have the same sizes as the element
   numbers in the mesh for element types of dimension ``spatial_dimension``.

If some changes have to be applied on the data as for example a padding for
``ParaView`` vectors, this can be done by using the field interface.

All these functions use the default dumper registered in the mesh but also have
the ``ToDumper`` variation with the dumper name specified. For example:

An example of code presenting this interface is present in the
``examples/io/dumper``. This interface is part of the  :cpp:class:`Dumpable
<akantu::Dumpable>` class from which the :cpp:class:`Mesh <akantu::Mesh>`
inherits.

Creating a new dumper
^^^^^^^^^^^^^^^^^^^^^

You can also create you own dumpers, ``Akantu`` uses a third-party library in
order to write the output files, ``IOHelper``. ``Akantu`` supports the
``ParaView`` format and a Text format defined by ``IOHelper``.

This two files format are handled by the classes :cpp:class:`DumperParaview
<akantu::DumperParaview>` and :cpp:class:`DumperText <akantu::DumperText>`.

In order to use them you can instantiate on of this object in your code. This
dumper have a simple interface. You can register a mesh :cpp:func:`registerMesh
<akantu::DumperIOHelper::registerMesh>`, :cpp:func:`registerFilteredMesh
<akantu::DumperIOHelper::registerFilteredMesh>` or a field,
:cpp:class:`registerField <akantu::DumperIOHelper::registerField>`.

An example of code presenting this low level interface is present in the
``examples/io/dumper``. The different types of :cpp:class:`Field
<akantu::Field>` that can be created are present in the source folder
``src/io/dumper``.
