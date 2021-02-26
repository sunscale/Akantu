Getting Started
===============

Compiling Akantu
----------------

Akantu is a `CMake <https://cmake.org/>`_ project, so to configure it, you can either
follow the usual way::

  > cd akantu
  > mkdir build
  > cd build
  > ccmake ..
  [ Set the options that you need ]
  > make
  > make install

All the Akantu options are documented in Appendix app:package-dependencies.

Writing a ``main`` function
---------------------------

Akantu first needs to be initialized. The memory management included in the core
library handles the correct allocation and de-allocation of vectors, structures
and/or objects. Moreover, in parallel computations, the initialization procedure
performs the communication setup. This is achieved by the function
:cpp:func:`initialize <akantu::initialize>` that is used as follows::

    #include "aka_common.hh"
    #include "..."

    using namespace akantu;

    int main(int argc, char *argv[]) {
	initialize("input_file.dat", argc, argv);

	// your code ...

    }

The :cpp:func:`initialize <akantu::initialize>` function takes the text inpute
file and the program parameters which can be parsed by Akantu in due form (see
sect:parser). Obviously it is necessary to include all files needed in main. In
this manual all provided code implies the usage of ``akantu`` as
namespace.

Creating and Loading a Mesh
---------------------------

In its current state, Akantu supports three types of meshes: Gmsh, Abaqus and
Diana. Once a :cpp:class:`akantu::Mesh` object is created with a given spatial
dimension, it can be filled by reading a mesh input file. The method
:cpp:func:`read <akantu::Mesh::read>` of the class :cpp:class:`Mesh
<akantu::Mesh>` infers the mesh type from the file extension. If a non-standard
file extension is used, the mesh type has to be specified. ::

    UInt spatial_dimension = 2;
    Mesh mesh(spatial_dimension);

    // Reading Gmsh files
    mesh.read("my_gmsh_mesh.msh");
    mesh.read("my_gmsh_mesh", _miot_gmsh);

The Gmsh reader adds the geometrical and physical tags as mesh data. The
physical values are stored as a :cpp:type:`UInt <akantu::UInt>` data called
``tag_0``, if a string name is provided it is stored as a ``std::string`` data
named ``physical_names``. The geometrical tag is stored as a :cpp:type:`UInt
<akantu::UInt>` data named ``tag_1``.

Using Arrays
------------

Data in Akantu can be stored in data containers implemented by the
:cpp:class:`akantu::Array` class. In its most basic usage, the :cpp:class:`Array
<akantu::Array>` class implemented in \akantu is similar to the ``std::vector``
class of the Standard Template Library (STL) for C++. A simple :cpp:class:`Array
<akantu::Array>` containing a sequence of ``nb_element`` values (of a given
type) can be generated with::

  Array<type> example_array(nb_element);

where ``type`` usually is ``Real``, ``Int``, ``UInt`` or ``bool``.
Each value is associated to an index, so that data can be accessed by
typing::

  auto & val = example_array(index);

``Arrays`` can also contain tuples of values for each index. In that case, the
number of components per tuple must be specified at the :cpp:class:`Array
<akantu::Array>` creation. For example, if we want to create an
:cpp:class:`Array <akantu::Array>` to store the coordinates (sequences of three
values) of ten nodes, the appropriate code is the following::

  UInt nb_nodes = 10;
  UInt spatial_dimension = 3;

  Array<Real> position(nb_nodes, spatial_dimension);

In this case the :math:`x` position of the eighth node number will be given
by ``position(7, 0)`` (in C++, numbering starts at 0 and not 1). If
the number of components for the sequences is not specified, the
default value of 1 is used. Here is a list of some basic operations
that can be performed on :cpp:class:`Array <akantu::Array>`:

  - ``resize(size)`` change the size of the :cpp:class:`Array <akantu::Array>`.
  - ``clear()`` set all entries of the :cpp:class:`Array <akantu::Array>` to
    zero.
  - ``set(t)`` set all entries of the :cpp:class:`Array <akantu::Array>` to
    ``t``.
  - ``copy(const Array<T> & other)`` copy another :cpp:class:`Array
    <akantu::Array>` into the current one. The two :cpp:class:`Array
    <akantu::Array>` should have the same number of components.
  - ``push_back(tuple)`` append a tuple with the correct number of components at
    the end of the :cpp:class:`Array <akantu::Array>`.
  - ``erase(i)`` erase the value at the i-th position.
  - ``find(value)`` search ``value`` in the current :cpp:class:`Array
    <akantu::Array>`. Return position index of the first occurence or -1 if not
    found.
  - ``storage()`` Return the address of the allocated memory of the
    :cpp:class:`Array <akantu::Array>`.

Array iterators
-------------------

It is very common in Akantu to loop over arrays to perform a specific treatment.
This ranges from geometric calculation on nodal quantities to tensor algebra (in
constitutive laws for example). The :cpp:class:`Array <akantu::Array>` object
has the possibility to request iterators in order to make the writing of loops
easier and enhance readability. For instance, a loop over the nodal coordinates
can be performed like::

  // accessing the nodal coordinates Array
  // with spatial_dimension components
  const auto & nodes = mesh.getNodes();

  for (const auto & coords : make_view(nodes, spatial_dimension)) {
    // do what you need ....
  }

In that example, each ``coords`` is a ``Vector<Real>`` containing
geometrical array of size ``spatial_dimension`` and the iteration is
conveniently performed by the :cpp:class:`Array <akantu::Array>` iterator.

The :cpp:class:`Array <akantu::Array>` object is intensively used to store
second order tensor values. In that case, it should be specified that the
returned object type is a matrix when constructing the iterator. This is done
when calling the :cpp:func:`make_view <akantu::make_view>`. For instance,
assuming that we have a :cpp:class:`Array <akantu::Array>` storing stresses, we
can loop over the stored tensors by::

   for (const auto & stress :
     make_view(stresses, spatial_dimension, spatial_dimension)) {
     // stress is of type `const Matrix<Real>&`
   }

In that last example, the :cpp:class:`Matrix <akantu::Matrix>` objects are
``spatial_dimension`` :math:`\times` ``spatial_dimension`` matrices. The light
objects :cpp:class:`Matrix <akantu::Matrix>` and :cpp:class:`Vector
<akantu::Vector>` can be used and combined to do most common linear algebra. If
the number of component is 1, it is possible to use :cpp:func:`make_view
<akantu::make_view>` to this effect.


In general, a mesh consists of several kinds of elements. Consequently, the
amount of data to be stored can differ for each element type. The
straightforward example is the connectivity array, namely the sequences of nodes
belonging to each element (linear triangular elements have fewer nodes than,
say, rectangular quadratic elements etc.). A particular data structure called
:cpp:class:`ElementTypeMapArray <akantu::ElementTypeMapArray>` is provided to
easily manage this kind of data. It consists of a group of ``Arrays``, each
associated to an element type. The following code can retrieve the
``ElementTypeMapArray`` which stores the connectivity arrays for a mesh::

  const ElementTypeMapArray<UInt> & connectivities =
    mesh.getConnectivities();

Then, the specific array associated to a given element type can be obtained by::

  const Array<UInt> & connectivity_triangle =
    connectivities(_triangle_3);

where the first order 3-node triangular element was used in the presented piece
of code.

Vector & Matrix
```````````````

The :cpp:class:`Array <akantu::Array>` iterators as presented in the previous
section can be shaped as :cpp:class:`Vector <akantu::Vector>` or
:cpp:class:`Matrix <akantu::Matrix>`. This objects represent 1st and 2nd order
tensors. As such they come with some functionalities that we will present a bit
more into detail in this here.


``Vector<T>``
'''''''''''''

- Accessors:

  - ``v(i)`` gives the ``i`` -th component of the vector ``v``
  - ``v[i]`` gives the ``i`` -th component of the vector ``v``
  - ``v.size()`` gives the number of component

- Level 1: (results are scalars)

  - ``v.norm()`` returns the geometrical norm (:math:`L_2`)
  - ``v.norm<N>()`` returns the :math:`L_N` norm defined as :math:`\left(\sum_i
    |v(i)|^N\right)^{1/N}`. N can take any positive integer value.
    There are also some particular values for the most commonly used
    norms, ``L_1`` for the Manhattan norm, ``L_2`` for the geometrical
    norm and ``L_inf`` for the norm infinity.
  - ``v.dot(x)`` return the dot product of ``v`` and ``x``
  - ``v.distance(x)`` return the geometrical norm of :math:`v - x`

- Level 2: (results are vectors)

  - ``v += s``, ``v -= s``, ``v *= s``, ``v /= s`` those are
    element-wise operators that sum, substract, multiply or divide all the
    component of ``v`` by the scalar ``s``
  - ``v += x``, ``v -= x`` sums or substracts the vector ``x`` to/from ``v``
  - ``v.mul(A, x, alpha)`` stores the result of :math:`\alpha \boldsymbol{A}
  \vec{x}` in ``v``, :math:`\alpha` is equal to 1 by default
  - ``v.solve(A, b)`` stores the result of the resolution of the system
    :math:`\boldsymbol{A} \vec{x} = \vec{b}` in ``v``
  - ``v.crossProduct(v1, v2)`` computes the cross product of ``v1`` and ``v2``
    and stores the result in ``v``

``Matrix<T>``
'''''''''''''

- Accessors:

  - ``A(i, j)`` gives the component :math:`A_{ij}` of the matrix ``A``
  - ``A(i)`` gives the :math:`i^{th}` column of the matrix as a ``Vector``
  - ``A[k]`` gives the :math:`k^{th}` component of the matrix, matrices are
    stored in a column major way, which means that to access :math:`A_{ij}`,
    :math:`k = i + j M`
  - ``A.rows()`` gives the number of rows of ``A`` (:math:`M`)
  - ``A.cols()`` gives the number of columns of ``A`` (:math:`N`)
  - ``A.size()`` gives the number of component in the matrix (:math:`M \times
    N`)

- Level 1: (results are scalars)

  - ``A.norm()`` is equivalent to ``A.norm<L_2>()``
  - ``A.norm<N>()`` returns the :math:`L_N` norm defined as
    :math:`\left(\sum_i\sum_j |A(i,j)|^N\right)^{1/N}`. N can take
    any positive integer value. There are also some particular values
    for the most commonly used norms, ``L_1`` for the Manhattan
    norm, ``L_2`` for the geometrical norm and ``L_inf`` for
    the norm infinity.
  - ``A.trace()`` return the trace of ``A``
  - ``A.det()`` return the determinant of ``A``
  - ``A.doubleDot(B)`` return the double dot product of ``A`` and
    ``B``, :math:`\mat{A}:\mat{B}`

- Level 3: (results are matrices)

  - ``A.eye(s)``, ``Matrix<T>::eye(s)`` fills/creates a matrix with
    the :math:`s\mat{I}` with :math:`\mat{I}` the identity matrix
  - ``A.inverse(B)`` stores :math:`\mat{B}^{-1}` in ``A``
  - ``A.transpose()`` returns  :math:`\mat{A}^{t}`
  - ``A.outerProduct(v1, v2)`` stores :math:`\vec{v_1} \vec{v_2}^{t}` in
    ``A``
  - ``C.mul<t_A, t_B>(A, B, alpha)``: stores the result of the product of
    ``A`` and code{B} time the scalar ``alpha`` in ``C``. ``t_A``
    and ``t_B`` are boolean defining if ``A`` and ``B`` should be
    transposed or not.

    +----------+----------+--------------+
    |``t_A``   |``t_B``   |result        |
    |          |          |              |
    +----------+----------+--------------+
    |false     |false     |:math:`\mat{C}|
    |          |          |= \alpha      |
    |          |          |\mat{A}       |
    |          |          |\mat{B}`      |
    |          |          |              |
    +----------+----------+--------------+
    |false     |true      |:math:`\mat{C}|
    |          |          |= \alpha      |
    |          |          |\mat{A}       |
    |          |          |\mat{B}^t`    |
    |          |          |              |
    +----------+----------+--------------+
    |true      |false     |:math:`\mat{C}|
    |          |          |= \alpha      |
    |          |          |\mat{A}^t     |
    |          |          |\mat{B}`      |
    |          |          |              |
    +----------+----------+--------------+
    |true      |true      |:math:`\mat{C}|
    |          |          |= \alpha      |
    |          |          |\mat{A}^t     |
    |          |          |\mat{B}^t`    |
    +----------+----------+--------------+

  - ``A.eigs(d, V)`` this method computes the eigenvalues and eigenvectors of
    ``A`` and store the results in ``d`` and ``V`` such that :math:`d(i) =
    \lambda_i` and :math:`V(i) = \vec{v_i}` with :math:`\mat{A}\vec{v_i} =
    \lambda_i\vec{v_i}` and :math:`\lambda_1 > ... > \lambda_i > ... >
    \lambda_N`
