``FEEngine``
============

The ``FEEngine`` interface is dedicated to handle the
finite-element approximations and the numerical integration of the
weak form. As we will see in Chapter sect:smm, ``Model``
creates its own ``FEEngine`` object so the explicit creation of the
object is not required.

Mathematical Operations
-----------------------

Using the ``FEEngine`` object, one can compute a interpolation, an
integration or a gradient. A simple example is given below::

    // having a FEEngine object
    std::unique_ptr<FEEngine> fem =
        std::make_unique<FEEngineTemplate<IntegratorGauss,ShapeLagrange>>(
            my_mesh, dim, "my_fem");
    // instead of this, a FEEngine object can be get using the model:
    // model.getFEEngine()

    //compute the gradient
    Array<Real> u; //append the values you want
    Array<Real> nablauq; //gradient array to be computed
    // compute the gradient
    fem->gradientOnIntegrationPoints(const Array<Real> &u,
                                     Array<Real> &nablauq,
                                     const UInt nb_degree_of_freedom,
                                     const ElementType & type);

    // interpolate
    Array<Real> uq; //interpolated array to be computed
    // compute the interpolation
    fem->interpolateOnIntegrationPoints(const Array<Real> &u,
                                        Array<Real> &uq,
                                        UInt nb_degree_of_freedom,
                                        const ElementType & type);

    // interpolated function can be integrated over the elements
    Array<Real> int_val_on_elem;
    // integrate
    fem->integrate(const Array<Real> &uq,
                   Array<Real> &int_uq,
                   UInt nb_degree_of_freedom,
                   const ElementType & type);

Another example below shows how to integrate stress and strain fields
over elements assigned to a particular material::

    UInt sp_dim = 3; //spatial dimension
    UInt m = 1; //material index of interest
    const ElementType type = _tetrahedron_4; //element type

    // get the stress and strain arrays associated to the material index m
    const Array<Real> & strain_vec = model.getMaterial(m).getGradU(type);
    const Array<Real> & stress_vec = model.getMaterial(m).getStress(type);

    // get the element filter for the material index
    const Array<UInt> & elem_filter = model.getMaterial(m).getElementFilter(type);

    // initialize the integrated stress and strain arrays
    Array<Real> int_strain_vec(elem_filter.getSize(),
                            sp_dim*sp_dim, "int_of_strain");
    Array<Real> int_stress_vec(elem_filter.getSize(),
                            sp_dim*sp_dim, "int_of_stress");

    // integrate the fields
    model.getFEEngine().integrate(strain_vec, int_strain_vec,
                                sp_dim*sp_dim, type, _not_ghost, elem_filter);
    model.getFEEngine().integrate(stress_vec, int_stress_vec,
                                sp_dim*sp_dim, type, _not_ghost, elem_filter);

Elements
--------

The base for every Finite-Elements computation is its mesh and the elements that
are used within that mesh. The element types that can be used depend on the
mesh, but also on the dimensionality of the problem (1D, 2D or 3D). In Akantu,
several isoparametric Lagrangian element types are supported (and one
serendipity element). Each of these types is discussed in some detail below,
starting with the 1D-elements all the way to the 3D-elements. More detailed
information (shape function, location of Gaussian quadrature points, and so on)
can be found in Appendix app:elements.

Isoparametric Elements
......................

1D
````

In Akantu, there are two types of isoparametric elements defined in 1D. These
element types are called ``_segment_2`` and ``_segment_3``, and are
depicted schematically in :numref:`fig:elements:1D`. Some of the basic
properties of these elements are listed in :numref:`tab:elements:1D`.

.. _fig:elements:1D:
.. figure:: figures/elements/segments.svg
            :align: center

            Schematic overview of the two 1D element types in Akantu. In each
            element, the node numbering as used in Akantu is indicated and also the
            quadrature points are highlighted (gray circles).


.. _tab:elements:1D:
.. table:: Some basic properties of the two 1D isoparametric elements in Akantu

   +--------------+---------+------+------+
   |Element       |Order    |#nodes|#quad |
   |type          |         |      |points|
   +--------------+---------+------+------+
   |``_segment_2``|linear   |2     |1     |
   +--------------+---------+------+------+
   |``_segment_3``|quadratic|3     |2     |
   +--------------+---------+------+------+

2D
````
In Akantu, there are four types of isoparametric elements defined in 2D. These
element types are called ``_triangle_3``, ``_triangle_6``,
``_quadrangle_4`` and ``_quadrangle_8``, and all of them are depicted
in :numref:`fig:elements:2D`. As with the 1D elements, some of the most basic
properties of these elements are listed in :numref:`tab:elements:2D`. It is
important to note that the first element is linear, the next two quadratic and
the last one cubic. Furthermore, the last element type (``_quadrangle_8``)
is not a Lagrangian but a serendipity element.

.. _fig:elements:2D:
.. figure:: figures/elements/elements_2d.svg
            :align: center

            Schematic overview of the four 2D element types in Akantu. In each
            element, the node numbering as used in Akantu is indicated and also the
            quadrature points are highlighted (gray circles).


.. _tab:elements:2D:
.. table:: Some basic properties of the 2D isoparametric elements in Akantu

   +--------------------+----------+------+------+
   |Element             |Order     |#nodes|#quad |
   |type                |          |      |points|
   +--------------------+----------+------+------+
   |``_triangle_3``     |linear    |3     |1     |
   +--------------------+----------+------+------+
   |``_triangle_6``     |quadratic |6     |3     |
   +--------------------+----------+------+------+
   |``_quadrangle_4``   |linear    |4     |4     |
   +--------------------+----------+------+------+
   |``_quadrangle_8``   |quadratic |8     |9     |
   +--------------------+----------+------+------+

3D
````

In Akantu, there are three types of isoparametric elements defined in 3D. These
element types are called ``_tetrahedron_4``, ``_tetrahedron_10`` and
``_hexahedron_8``, and all of them are depicted schematically in
:numref:`fig:elements:3D`. As with the 1D and 2D elements some of the most
basic properties of these elements are listed in :numref:`tab:elements:3D`.

.. _fig:elements:3D:
.. figure:: figures/elements/elements_3d.svg
            :align: center

            Schematic overview of the three 3D element types in Akantu. In each
            element, the node numbering as used in Akantu is indicated and also the
            quadrature points are highlighted (gray circles).

.. _tab:elements:3D:
.. table:: Some basic properties of the 3D isoparametric elements in Akantu

   +--------------------+----------+------+------+
   |Element             |Order     |#nodes|#quad |
   |type                |          |      |points|
   +--------------------+----------+------+------+
   |``_tetrahedron_4``  |linear    |4     |1     |
   +--------------------+----------+------+------+
   |``_tetrahedron_10`` |quadratic |10    |4     |
   +--------------------+----------+------+------+
   |``_hexadedron_8``   |cubic     |8     |8     |
   +--------------------+----------+------+------+
