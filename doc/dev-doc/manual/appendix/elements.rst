.. _app-elements:

Shape Functions
===============

Schematic overview of all the element types defined in `Akantu` is described in
Section :ref:`sec-elements`. In this appendix, more detailed information (shape
function, location of Gaussian quadrature points, and so on) of each of these
types is listed. For each element type, the coordinates of the nodes are given
in the iso-parametric frame of reference, together with the shape functions (and
their derivatives) on these respective nodes. Also all the Gaussian quadrature
points within each element are assigned (together with the weight that is
applied on these points). The graphical representations of all the element types
can be found in Section :ref:`sec-elements`.

Iso-parametric Elements
-----------------------

1D-Shape Functions
``````````````````

Segment 2
'''''''''

.. list-table:: Elements properties
   :header-rows: 1

   * - Node (:math:`i`)
     - Coord. (:math:`\xi`)
     - Shape function (:math:`N_i`)
     - Derivative (:math:`\frac{\partial N_i}{\partial \xi}`)
   * - 1
     - -1
     - :math:`\frac{1}{2}\left(1-\xi\right)`
     - :math:`-\frac{1}{2}`
   * - 2
     - 1
     - :math:`\frac{1}{2}\left(1+\xi\right)`
     - :math:`\frac{1}{2}`

.. list-table:: Gaussian quadrature points
   :align: center

   * - Coord. (:math:`\xi`)
     - Weight
   * - 0
     - 2

Segment 3
'''''''''

.. list-table:: Elements properties
   :header-rows: 1

   * - Node (:math:`i`)
     - Coord. (:math:`\xi`)
     - Shape function (:math:`N_i`)
     - Derivative (:math:`\frac{\partial N_i}{\partial \xi}`)
   * - 1
     - -1
     - :math:`\frac{1}{2}\xi\left(\xi-1\right)`
     - :math:`\xi-\frac{1}{2}`
   * - 2
     - 1
     - :math:`\frac{1}{2}\xi\left(\xi+1\right)`
     - :math:`\xi+\frac{1}{2}`
   * - 3
     - 0
     - :math:`1-\xi^{2}`
     - :math:`-2\xi`

.. list-table:: Gaussian quadrature points
   :align: center

   * - Coord. (:math:`\xi`)
     - Weight
   * - :math:`-1/\sqrt{3}`
     - 1
   * - :math:`1/\sqrt{3}`
     - 1


2D-Shape Functions
``````````````````

Triangle 3
''''''''''

.. list-table:: Elements properties
   :header-rows: 1

   * - Node (:math:`i`)
     - Coord. (:math:`\xi`, :math:`\eta`)
     - Shape function (:math:`N_i`)
     - Derivative (:math:`\frac{\partial N_i}{\partial \xi}`, :math:`\frac{\partial N_i}{\partial \eta}`)
   * - 1
     - (:math:`0`, :math:`0`)
     - :math:`1-\xi-\eta`
     - (:math:`-1`, :math:`-1`)
   * - 2
     - (:math:`1`, :math:`0`)
     - :math:`\xi`
     - (:math:`1`, :math:`0`)
   * - 3
     - (:math:`0`, :math:`1`)
     - :math:`\eta`
     - (:math:`0`, :math:`1`)

.. list-table:: Gaussian quadrature points
   :align: center

   * - Coord. (:math:`\xi`, :math:`\eta`)
     - Weight
   * - (:math:`\frac{1}{3}`, :math:`\frac{1}{3}`)
     - :math:`\frac{1}{2}`

Triangle 6
''''''''''

.. list-table:: Elements properties
   :header-rows: 1

   * - Node (:math:`i`)
     - Coord. (:math:`\xi`, :math:`\eta`)
     - Shape function (:math:`N_i`)
     - Derivative (:math:`\frac{\partial N_i}{\partial \xi}`, :math:`\frac{\partial N_i}{\partial \eta}`)
   * - 1
     - (:math:`0`, :math:`0`)
     - :math:`-\left(1-\xi-\eta\right)\left(1-2\left(1-\xi-\eta\right)\right)`
     - (:math:`1-4\left(1-\xi-\eta\right)`, :math:`1-4\left(1-\xi-\eta\right)`)
   * - 2
     - (:math:`1`, :math:`0`)
     - :math:`-\xi\left(1-2\xi\right)`
     - (:math:`4\xi-1`, :math:`0`)
   * - 3
     - (:math:`0`, :math:`1`)
     - :math:`-\eta\left(1-2\eta\right)`
     - (:math:`0`, :math:`4\eta-1`)
   * - 4
     - (:math:`\frac{1}{2}`, :math:`0`)
     - :math:`4\xi\left(1-\xi-\eta\right)`
     - (:math:`4\left(1-2\xi-\eta\right)`, :math:`-4\xi`)
   * - 5
     - (:math:`\frac{1}{2}`, :math:`\frac{1}{2}`)
     - :math:`4\xi\eta`
     - (:math:`4\eta`, :math:`4\xi`)
   * - 6
     - (:math:`0`, :math:`\frac{1}{2}`)
     - :math:`4\eta\left(1-\xi-\eta\right)`
     - (:math:`-4\eta`, :math:`4\left(1-\xi-2\eta\right)`)

.. list-table:: Gaussian quadrature points
   :align: center


   * - Coord. (:math:`\xi`, :math:`\eta`)
     - Weight
   * - (:math:`\frac{1}{6}`, :math:`\frac{1}{6}`)
     - :math:`\frac{1}{6}`
   * - (:math:`\frac{2}{3}`, :math:`\frac{1}{6}`)
     - :math:`\frac{1}{6}`
   * - (:math:`\frac{1}{6}`, :math:`\frac{2}{3}`)
     - :math:`\frac{1}{6}`

Quadrangle 4
''''''''''''

.. list-table:: Elements properties
   :header-rows: 1

   * - Node (:math:`i`)
     - Coord. (:math:`\xi`, :math:`\eta`)
     - Shape function (:math:`N_i`)
     - Derivative (:math:`\frac{\partial N_i}{\partial \xi}`, :math:`\frac{\partial N_i}{\partial \eta}`)
   * - 1
     - (:math:`-1`, :math:`-1`)
     - :math:`\frac{1}{4}\left(1-\xi\right)\left(1-\eta\right)`
     - (:math:`-\frac{1}{4}\left(1-\eta\right)`, :math:`-\frac{1}{4}\left(1-\xi\right)`)
   * - 2
     - (:math:`1`, :math:`-1`)
     - :math:`\frac{1}{4}\left(1+\xi\right)\left(1-\eta\right)`
     - (:math:`\frac{1}{4}\left(1-\eta\right)`, :math:`-\frac{1}{4}\left(1+\xi\right)`)
   * - 3
     - (:math:`1`, :math:`1`)
     - :math:`\frac{1}{4}\left(1+\xi\right)\left(1+\eta\right)`
     - (:math:`\frac{1}{4}\left(1+\eta\right)`, :math:`\frac{1}{4}\left(1+\xi\right)`)
   * - 4
     - (:math:`-1`, :math:`1`)
     - :math:`\frac{1}{4}\left(1-\xi\right)\left(1+\eta\right)`
     - (:math:`-\frac{1}{4}\left(1+\eta\right)`, :math:`\frac{1}{4}\left(1-\xi\right)`)

.. list-table:: Gaussian quadrature points
   :align: center

   * - Coord. (:math:`\xi`, :math:`\eta`)
     - Weight
   * - (:math:`-\frac{1}{\sqrt{3}}`, :math:`-\frac{1}{\sqrt{3}}`)
     - 1
   * - (:math:`\frac{1}{\sqrt{3}}`, :math:`-\frac{1}{\sqrt{3}}`)
     - 1
   * - (:math:`\frac{1}{\sqrt{3}}`, :math:`\frac{1}{\sqrt{3}}`)
     - 1
   * - (:math:`-\frac{1}{\sqrt{3}}`, :math:`\frac{1}{\sqrt{3}}`)
     - 1

Quadrangle 8
''''''''''''

.. list-table:: Elements properties
   :header-rows: 1

   * - Node (:math:`i`)
     - Coord. (:math:`\xi`, :math:`\eta`)
     - Shape function (:math:`N_i`)
     - Derivative (:math:`\frac{\partial N_i}{\partial \xi}`, :math:`\frac{\partial N_i}{\partial \eta}`)
   * - 1
     - (:math:`-1`, :math:`-1`)
     - :math:`\frac{1}{4}\left(1-\xi\right)\left(1-\eta\right)\left(-1-\xi-\eta\right)`
     - (:math:`\frac{1}{4}\left(1-\eta\right)\left(2\xi+\eta\right)`, :math:`\frac{1}{4}\left(1-\xi\right)\left(\xi+2\eta\right)`)
   * - 2
     - (:math:`1`, :math:`-1`)
     - :math:`\frac{1}{4}\left(1+\xi\right)\left(1-\eta\right)\left(-1+\xi-\eta\right)`
     - (:math:`\frac{1}{4}\left(1-\eta\right)\left(2\xi-\eta\right)`, :math:`-\frac{1}{4}\left(1+\xi\right)\left(\xi-2\eta\right)`)
   * - 3
     - (:math:`1`, :math:`1`)
     - :math:`\frac{1}{4}\left(1+\xi\right)\left(1+\eta\right)\left(-1+\xi+\eta\right)`
     - (:math:`\frac{1}{4}\left(1+\eta\right)\left(2\xi+\eta\right)`, :math:`\frac{1}{4}\left(1+\xi\right)\left(\xi+2\eta\right)`)
   * - 4
     - (:math:`-1`, :math:`1`)
     - :math:`\frac{1}{4}\left(1-\xi\right)\left(1+\eta\right)\left(-1-\xi+\eta\right)`
     - (:math:`\frac{1}{4}\left(1+\eta\right)\left(2\xi-\eta\right)`, :math:`-\frac{1}{4}\left(1-\xi\right)\left(\xi-2\eta\right)`)
   * - 5
     - (:math:`0`, :math:`-1`)
     - :math:`\frac{1}{2}\left(1-\xi^{2}\right)\left(1-\eta\right)`
     - (:math:`-\xi\left(1-\eta\right)`, :math:`-\frac{1}{2}\left(1-\xi^{2}\right)`)
   * - 6
     - (:math:`1`, :math:`0`)
     - :math:`\frac{1}{2}\left(1+\xi\right)\left(1-\eta^{2}\right)`
     - (:math:`\frac{1}{2}\left(1-\eta^{2}\right)`, :math:`-\eta\left(1+\xi\right)`)
   * - 7
     - (:math:`0`, :math:`1`)
     - :math:`\frac{1}{2}\left(1-\xi^{2}\right)\left(1+\eta\right)`
     - (:math:`-\xi\left(1+\eta\right)`, :math:`\frac{1}{2}\left(1-\xi^{2}\right)`)
   * - 8
     - (:math:`-1`, :math:`0`)
     - :math:`\frac{1}{2}\left(1-\xi\right)\left(1-\eta^{2}\right)`
     - (:math:`-\frac{1}{2}\left(1-\eta^{2}\right)`, :math:`-\eta\left(1-\xi\right)`)

.. list-table:: Gaussian quadrature points
   :align: center

   * - Coord. (:math:`\xi`, :math:`\eta`)
     - Weight
   * - (:math:`0`, :math:`0`)
     - :math:`\frac{64}{81}`
   * - (:math:`\sqrt{\tfrac{3}{5}}`, :math:`\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{25}{81}`
   * - (:math:`-\sqrt{\tfrac{3}{5}}`, :math:`\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{25}{81}`
   * - (:math:`-\sqrt{\tfrac{3}{5}}`, :math:`-\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{25}{81}`
   * - (:math:`\sqrt{\tfrac{3}{5}}`, :math:`-\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{25}{81}`
   * - (:math:`0`, :math:`\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{40}{81}`
   * - (:math:`-\sqrt{\tfrac{3}{5}}`, :math:`0`)
     - :math:`\frac{40}{81}`
   * - (:math:`0`, :math:`-\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{40}{81}`
   * - (:math:`\sqrt{\tfrac{3}{5}}`, :math:`0`)
     - :math:`\frac{40}{81}`


3D-Shape Functions
``````````````````

Tetrahedron 4
'''''''''''''

.. list-table:: Elements properties
   :header-rows: 1

   * - Node (:math:`i`)
     - Coord. (:math:`\xi`, :math:`\eta`, :math:`\zeta`)
     - Shape function (:math:`N_i`)
     - Derivative (:math:`\frac{\partial N_i}{\partial \xi}`, :math:`\frac{\partial N_i}{\partial \eta}`, :math:`\frac{\partial N_i}{\partial \zeta}`)
   * - 1
     - (:math:`0`, :math:`0`, :math:`0`)
     - :math:`1-\xi-\eta-\zeta`
     - (:math:`-1`, :math:`-1`, :math:`-1`)
   * - 2
     - (:math:`1`, :math:`0`, :math:`0`)
     - :math:`\xi`
     - (:math:`1`, :math:`0`, :math:`0`)
   * - 3
     - (:math:`0`, :math:`1`, :math:`0`)
     - :math:`\eta`
     - (:math:`0`, :math:`1`, :math:`0`)
   * - 4
     - (:math:`0`, :math:`0`, :math:`1`)
     - :math:`\zeta`
     - (:math:`0`, :math:`0`, :math:`1`)

.. list-table:: Gaussian quadrature points
   :align: center

   * - Coord. (:math:`\xi`, :math:`\eta`, :math:`\zeta`)
     - Weight
   * - (:math:`\frac{1}{4}`, :math:`\frac{1}{4}`, :math:`\frac{1}{4}`)
     - :math:`\frac{1}{6}`

Tetrahedron 10
''''''''''''''

.. list-table:: Elements properties
   :header-rows: 1

   * - Node (:math:`i`)
     - Coord. (:math:`\xi`, :math:`\eta`, :math:`\zeta`)
     - Shape function (:math:`N_i`)
     - Derivative (:math:`\frac{\partial N_i}{\partial \xi}`, :math:`\frac{\partial N_i}{\partial \eta}`, :math:`\frac{\partial N_i}{\partial \zeta}`)
   * - 1
     - (:math:`0`, :math:`0`, :math:`0`)
     - :math:`\left(1-\xi-\eta-\zeta\right)\left(1-2\xi-2\eta-2\zeta\right)`
     - :math:`4\xi+4\eta+4\zeta-3`, :math:`4\xi+4\eta+4\zeta-3`, :math:`4\xi+4\eta+4\zeta-3`
   * - 2
     - (:math:`1`, :math:`0`, :math:`0`)
     - :math:`\xi\left(2\xi-1\right)`
     - (:math:`4\xi-1`, :math:`0`, :math:`0`)
   * - 3
     - (:math:`0`, :math:`1`, :math:`0`)
     - :math:`\eta\left(2\eta-1\right)`
     - (:math:`0`, :math:`4\eta-1`, :math:`0`)
   * - 4
     - (:math:`0`, :math:`0`, :math:`1`)
     - :math:`\zeta\left(2\zeta-1\right)`
     - (:math:`0`, :math:`0`, :math:`4\zeta-1`)
   * - 5
     - (:math:`\frac{1}{2}`, :math:`0`, :math:`0`)
     - :math:`4\xi\left(1-\xi-\eta-\zeta\right)`
     - (:math:`4-8\xi-4\eta-4\zeta`, :math:`-4\xi`, :math:`-4\xi`)
   * - 6
     - (:math:`\frac{1}{2}`, :math:`\frac{1}{2}`, :math:`0`)
     - :math:`4\xi\eta`
     - (:math:`4\eta`, :math:`4\xi`, :math:`0`)
   * - 7
     - (:math:`0`, :math:`\frac{1}{2}`, :math:`0`)
     - :math:`4\eta\left(1-\xi-\eta-\zeta\right)`
     - (:math:`-4\eta`, :math:`4-4\xi-8\eta-4\zeta`, :math:`-4\eta`)
   * - 8
     - (:math:`0`, :math:`0`, :math:`\frac{1}{2}`)
     - :math:`4\zeta\left(1-\xi-\eta-\zeta\right)`
     - (:math:`-4\zeta`, :math:`-4\zeta`, :math:`4-4\xi-4\eta-8\zeta`)
   * - 9
     - (:math:`\frac{1}{2}`, :math:`0`, :math:`\frac{1}{2}`)
     - :math:`4\xi\zeta`
     - (:math:`4\zeta`, :math:`0`, :math:`4\xi`)
   * - 10
     - (:math:`0`, :math:`\frac{1}{2}`, :math:`\frac{1}{2}`)
     - :math:`4\eta\zeta`
     - (:math:`0`, :math:`4\zeta`, :math:`4\eta`)

.. list-table:: Gaussian quadrature points
   :align: center

   * - Coord. (:math:`\xi`, :math:`\eta`, :math:`\zeta`)
     - Weight
   * - (:math:`\frac{5-\sqrt{5}}{20}`, :math:`\frac{5-\sqrt{5}}{20}`, :math:`\frac{5-\sqrt{5}}{20}`)
     - :math:`\frac{1}{24}`
   * - (:math:`\frac{5+3\sqrt{5}}{20}`, :math:`\frac{5-\sqrt{5}}{20}`, :math:`\frac{5-\sqrt{5}}{20}`)
     - :math:`\frac{1}{24}`
   * - (:math:`\frac{5-\sqrt{5}}{20}`, :math:`\frac{5+3\sqrt{5}}{20}`, :math:`\frac{5-\sqrt{5}}{20}`)
     - :math:`\frac{1}{24}`
   * - (:math:`\frac{5-\sqrt{5}}{20}`, :math:`\frac{5-\sqrt{5}}{20}`, :math:`\frac{5+3\sqrt{5}}{20}`)
     - :math:`\frac{1}{24}`

Hexahedron 8
''''''''''''

.. list-table:: Elements properties
   :header-rows: 1

   * - Node (:math:`i`)
     - Coord. (:math:`\xi`, :math:`\eta`, :math:`\zeta`)
     - Shape function (:math:`N_i`)
     - Derivative (:math:`\frac{\partial N_i}{\partial \xi}`, :math:`\frac{\partial N_i}{\partial \eta}`, :math:`\frac{\partial N_i}{\partial \zeta}`)
   * - 1
     - (:math:`-1`, :math:`-1`, :math:`-1`)
     - :math:`\frac{1}{8}\left(1-\xi\right)\left(1-\eta\right)\left(1-\zeta\right)`
     - (:math:`-\frac{1}{8}\left(1-\eta\right)\left(1-\zeta\right)`, :math:`-\frac{1}{8}\left(1-\xi\right)\left(1-\zeta\right)`, :math:`3`)
   * - 2
     - (:math:`1`, :math:`-1`, :math:`-1`)
     - :math:`\frac{1}{8}\left(1+\xi\right)\left(1-\eta\right)\left(1-\zeta\right)`
     - (:math:`\frac{1}{8}\left(1-\eta\right)\left(1-\zeta\right)`, :math:`-\frac{1}{8}\left(1+\xi\right)\left(1-\zeta\right)`, :math:`3`)
   * - 3
     - (:math:`1`, :math:`1`, :math:`-1`)
     - :math:`\frac{1}{8}\left(1+\xi\right)\left(1+\eta\right)\left(1-\zeta\right)`
     - (:math:`\frac{1}{8}\left(1+\eta\right)\left(1-\zeta\right)`, :math:`\frac{1}{8}\left(1+\xi\right)\left(1-\zeta\right)`, :math:`3`)
   * - 4
     - (:math:`-1`, :math:`1`, :math:`-1`)
     - :math:`\frac{1}{8}\left(1-\xi\right)\left(1+\eta\right)\left(1-\zeta\right)`
     - (:math:`-\frac{1}{8}\left(1+\eta\right)\left(1-\zeta\right)`, :math:`\frac{1}{8}\left(1-\xi\right)\left(1-\zeta\right)`, :math:`3`)
   * - 5
     - (:math:`-1`, :math:`-1`, :math:`1`)
     - :math:`\frac{1}{8}\left(1-\xi\right)\left(1-\eta\right)\left(1+\zeta\right)`
     - (:math:`-\frac{1}{8}\left(1-\eta\right)\left(1+\zeta\right)`, :math:`-\frac{1}{8}\left(1-\xi\right)\left(1+\zeta\right)`, :math:`3`)
   * - 6
     - (:math:`1`, :math:`-1`, :math:`1`)
     - :math:`\frac{1}{8}\left(1+\xi\right)\left(1-\eta\right)\left(1+\zeta\right)`
     - (:math:`\frac{1}{8}\left(1-\eta\right)\left(1+\zeta\right)`, :math:`-\frac{1}{8}\left(1+\xi\right)\left(1+\zeta\right)`, :math:`3`)
   * - 7
     - (:math:`1`, :math:`1`, :math:`1`)
     - :math:`\frac{1}{8}\left(1+\xi\right)\left(1+\eta\right)\left(1+\zeta\right)`
     - (:math:`\frac{1}{8}\left(1+\eta\right)\left(1+\zeta\right)`, :math:`\frac{1}{8}\left(1+\xi\right)\left(1+\zeta\right)`, :math:`3`)
   * - 8
     - (:math:`-1`, :math:`1`, :math:`1`)
     - :math:`\frac{1}{8}\left(1-\xi\right)\left(1+\eta\right)\left(1+\zeta\right)`
     - (:math:`-\frac{1}{8}\left(1+\eta\right)\left(1+\zeta\right)`, :math:`\frac{1}{8}\left(1-\xi\right)\left(1+\zeta\right)`, :math:`3`)

.. list-table:: Gaussian quadrature points
   :align: center

   * - Coord. (:math:`\xi`, :math:`\eta`, :math:`\zeta`)
     - Weight
   * - (:math:`-\frac{1}{\sqrt{3}}`, :math:`-\frac{1}{\sqrt{3}}`, :math:`-\frac{1}{\sqrt{3}}`)
     - 1
   * - (:math:`\frac{1}{\sqrt{3}}`, :math:`-\frac{1}{\sqrt{3}}`, :math:`-\frac{1}{\sqrt{3}}`)
     - 1
   * - (:math:`\frac{1}{\sqrt{3}}`, :math:`\frac{1}{\sqrt{3}}`, :math:`-\frac{1}{\sqrt{3}}`)
     - 1
   * - (:math:`-\frac{1}{\sqrt{3}}`, :math:`\frac{1}{\sqrt{3}}`, :math:`-\frac{1}{\sqrt{3}}`)
     - 1
   * - (:math:`-\frac{1}{\sqrt{3}}`, :math:`-\frac{1}{\sqrt{3}}`, :math:`\frac{1}{\sqrt{3}}`)
     - 1
   * - (:math:`\frac{1}{\sqrt{3}}`, :math:`-\frac{1}{\sqrt{3}}`, :math:`\frac{1}{\sqrt{3}}`)
     - 1
   * - (:math:`\frac{1}{\sqrt{3}}`, :math:`\frac{1}{\sqrt{3}}`, :math:`\frac{1}{\sqrt{3}}`)
     - 1
   * - (:math:`-\frac{1}{\sqrt{3}}`, :math:`\frac{1}{\sqrt{3}}`, :math:`\frac{1}{\sqrt{3}}`)
     - 1

Pentahedron 6
'''''''''''''

.. list-table:: Elements properties
   :header-rows: 1

   * - Node (:math:`i`)
     - Coord. (:math:`\xi`, :math:`\eta`, :math:`\zeta`)
     - Shape function (:math:`N_i`)
     - Derivative (:math:`\frac{\partial N_i}{\partial \xi}`, :math:`\frac{\partial N_i}{\partial \eta}`, :math:`\frac{\partial N_i}{\partial \zeta}`)
   * - 1
     - (:math:`-1`, :math:`1`, :math:`0`)
     - :math:`\frac{1}{2}\left(1-\xi\right)\eta`

     - (:math:`-\frac{1}{2}\eta`, :math:`\frac{1}{2}\left(1-\xi\right)`, :math:`3`)
   * - 2
     - (:math:`-1`, :math:`0`, :math:`1`)
     - :math:`\frac{1}{2}\left(1-\xi\right)\zeta`

     - (:math:`-\frac{1}{2}\zeta`, :math:`0.0`, :math:`3`)
   * - 3
     - (:math:`-1`, :math:`0`, :math:`0`)
     - :math:`\frac{1}{2}\left(1-\xi\right)\left(1-\eta-\zeta\right)`

     - (:math:`-\frac{1}{2}\left(1-\eta-\zeta\right)`, :math:`-\frac{1}{2}\left(1-\xi\right)`, :math:`3`)
   * - 4
     - (:math:`1`, :math:`1`, :math:`0`)
     - :math:`\frac{1}{2}\left(1+\xi\right)\eta`

     - (:math:`\frac{1}{2}\eta`, :math:`\frac{1}{2}\left(1+\xi\right)`, :math:`3`)
   * - 5
     - (:math:`1`, :math:`0`, :math:`1`)
     - :math:`\frac{1}{2}\left(1+\xi\right)\zeta`

     - (:math:`\frac{1}{2}\zeta`, :math:`0.0`, :math:`3`)
   * - 6
     - (:math:`1`, :math:`0`, :math:`0`)
     - :math:`\frac{1}{2}\left(1+\xi\right)\left(1-\eta-\zeta\right)`

     - (:math:`\frac{1}{2}\left(1-\eta-\zeta\right)`, :math:`-\frac{1}{2}\left(1+\xi\right)`, :math:`3`)

.. list-table:: Gaussian quadrature points
   :align: center

   * - Coord. (:math:`\xi`, :math:`\eta`, :math:`\zeta`)
     - Weight
   * - (:math:`-\frac{1}{\sqrt{3}}`, :math:`0.5`, :math:`0.5`)
     - :math:`\frac{1}{6}`
   * - (:math:`-\frac{1}{\sqrt{3}}`, :math:`0.0`, :math:`0.5`)
     - :math:`\frac{1}{6}`
   * - (:math:`-\frac{1}{\sqrt{3}}`, :math:`0.5`, :math:`0.0`)
     - :math:`\frac{1}{6}`
   * - (:math:`\frac{1}{\sqrt{3}}`, :math:`0.5`, :math:`0.5`)
     - :math:`\frac{1}{6}`
   * - (:math:`\frac{1}{\sqrt{3}}`, :math:`0.0`, :math:`0.5`)
     - :math:`\frac{1}{6}`
   * - (:math:`\frac{1}{\sqrt{3}}`, :math:`0.5`, :math:`0.0`)
     - :math:`\frac{1}{6}`

Hexahedron 20
'''''''''''''

.. list-table:: Elements properties
   :header-rows: 1

   * - Node (:math:`i`)
     - Coord. (:math:`\xi`, :math:`\eta`, :math:`\zeta`)
     - Shape function (:math:`N_i`)
     - Derivative (:math:`\frac{\partial N_i}{\partial \xi}`, :math:`\frac{\partial N_i}{\partial \eta}`, :math:`\frac{\partial N_i}{\partial \zeta}`)
   * - 1
     - (:math:`-1`, :math:`-1`, :math:`-1`)
     - :math:`\frac{1}{8}\left(1-\xi\right)\left(1-\eta\right)\left(1-\zeta\right)\left(-2-\xi-\eta-\zeta\right)`
     - (:math:`\frac{1}{4}\left(\xi+\frac{1}{2}\left(\eta+\zeta+1\right)\right)\left(\eta-1\right)\left(\zeta-1\right)`, :math:`\frac{1}{4}\left(\eta+\frac{1}{2}\left(\xi+\zeta+1\right)\right)\left(\xi-1\right)\left(\zeta-1\right)`, :math:`3`)
   * - 2
     - (:math:`1`, :math:`-1`, :math:`-1`)
     - :math:`\frac{1}{8}\left(1+\xi\right)\left(1-\eta\right)\left(1-\zeta\right)\left(-2+\xi-\eta-\zeta\right)`
     - (:math:`\frac{1}{4}\left(\xi-\frac{1}{2}\left(\eta+\zeta+1\right)\right)\left(\eta-1\right)\left(\zeta-1\right)`, :math:`-\frac{1}{4}\left(\eta-\frac{1}{2}\left(\xi-\zeta-1\right)\right)\left(\xi+1\right)\left(\zeta-1\right)`, :math:`3`)
   * - 3
     - (:math:`1`, :math:`1`, :math:`-1`)
     - :math:`\frac{1}{8}\left(1+\xi\right)\left(1+\eta\right)\left(1-\zeta\right)\left(-2+\xi+\eta-\zeta\right)`
     - (:math:`-\frac{1}{4}\left(\xi+\frac{1}{2}\left(\eta-\zeta-1\right)\right)\left(\eta+1\right)\left(\zeta-1\right)`, :math:`-\frac{1}{4}\left(\eta+\frac{1}{2}\left(\xi-\zeta-1\right)\right)\left(\xi+1\right)\left(\zeta-1\right)`, :math:`3`)
   * - 4
     - (:math:`-1`, :math:`1`, :math:`-1`)
     - :math:`\frac{1}{8}\left(1-\xi\right)\left(1+\eta\right)\left(1-\zeta\right)\left(-2-\xi+\eta-\zeta\right)`
     - (:math:`-\frac{1}{4}\left(\xi-\frac{1}{2}\left(\eta-\zeta-1\right)\right)\left(\eta+1\right)\left(\zeta-1\right)`, :math:`\frac{1}{4}\left(\eta-\frac{1}{2}\left(\xi+\zeta+1\right)\right)\left(\xi-1\right)\left(\zeta-1\right)`, :math:`3`)
   * - 5
     - (:math:`-1`, :math:`-1`, :math:`1`)
     - :math:`\frac{1}{8}\left(1-\xi\right)\left(1-\eta\right)\left(1+\zeta\right)\left(-2-\xi-\eta+\zeta\right)`
     - (:math:`-\frac{1}{4}\left(\xi+\frac{1}{2}\left(\eta-\zeta+1\right)\right)\left(\eta-1\right)\left(\zeta+1\right)`, :math:`-\frac{1}{4}\left(\eta+\frac{1}{2}\left(\xi-\zeta+1\right)\right)\left(\xi-1\right)\left(\zeta+1\right)`, :math:`3`)
   * - 6
     - (:math:`1`, :math:`-1`, :math:`1`)
     - :math:`\frac{1}{8}\left(1+\xi\right)\left(1-\eta\right)\left(1+\zeta\right)\left(-2+\xi-\eta+\zeta\right)`
     - (:math:`-\frac{1}{4}\left(\xi-\frac{1}{2}\left(\eta-\zeta+1\right)\right)\left(\eta-1\right)\left(\zeta+1\right)`, :math:`\frac{1}{4}\left(\eta-\frac{1}{2}\left(\xi+\zeta-1\right)\right)\left(\xi+1\right)\left(\zeta+1\right)`, :math:`3`)
   * - 7
     - (:math:`1`, :math:`1`, :math:`1`)
     - :math:`\frac{1}{8}\left(1+\xi\right)\left(1+\eta\right)\left(1+\zeta\right)\left(-2+\xi+\eta+\zeta\right)`
     - (:math:`\frac{1}{4}\left(\xi+\frac{1}{2}\left(\eta+\zeta-1\right)\right)\left(\eta+1\right)\left(\zeta+1\right)`, :math:`\frac{1}{4}\left(\eta+\frac{1}{2}\left(\xi+\zeta-1\right)\right)\left(\xi+1\right)\left(\zeta+1\right)`, :math:`3`)
   * - 8
     - (:math:`-1`, :math:`1`, :math:`1`)
     - :math:`\frac{1}{8}\left(1-\xi\right)\left(1+\eta\right)\left(1+\zeta\right)\left(-2-\xi+\eta+\zeta\right)`
     - (:math:`\frac{1}{4}\left(\xi-\frac{1}{2}\left(\eta+\zeta-1\right)\right)\left(\eta+1\right)\left(\zeta+1\right)`, :math:`-\frac{1}{4}\left(\eta-\frac{1}{2}\left(\xi-\zeta+1\right)\right)\left(\xi-1\right)\left(\zeta+1\right)`, :math:`3`)
   * - 9
     - (:math:`0`, :math:`-1`, :math:`-1`)
     - :math:`\frac{1}{4}\left(1-\xi^{2}\right)\left(1-\eta\right)\left(1-\zeta\right)`
     - (:math:`-\frac{1}{2}\xi\left(\eta-1\right)\left(\zeta-1\right)`, :math:`-\frac{1}{4}\left(\xi^{2}-1\right)\left(\zeta-1\right)`, :math:`3`)
   * - 10
     - (:math:`1`, :math:`0`, :math:`-1`)
     - :math:`\frac{1}{4}\left(1+\xi\right)\left(1-\eta^{2}\right)\left(1-\zeta\right)`
     - (:math:`\frac{1}{4}\left(\eta^{2}-1\right)\left(\zeta-1\right)`, :math:`\frac{1}{2}\eta\left(\xi+1\right)\left(\zeta-1\right)`, :math:`3`)
   * - 11
     - (:math:`0`, :math:`1`, :math:`-1`)
     - :math:`\frac{1}{4}\left(1-\xi^{2}\right)\left(1+\eta\right)\left(1-\zeta\right)`
     - (:math:`\frac{1}{2}\xi\left(\eta+1\right)\left(\zeta-1\right)`, :math:`\frac{1}{4}\left(\xi^{2}-1\right)\left(\zeta-1\right)`, :math:`3`)
   * - 12
     - (:math:`-1`, :math:`0`, :math:`-1`)
     - :math:`\frac{1}{4}\left(1-\xi\right)\left(1-\eta^{2}\right)\left(1-\zeta\right)`
     - (:math:`-\frac{1}{4}\left(\eta^{2}-1\right)\left(\zeta-1\right)`, :math:`-\frac{1}{2}\eta\left(\xi-1\right)\left(\zeta-1\right)`, :math:`3`)
   * - 13
     - (:math:`-1`, :math:`-1`, :math:`0`)
     - :math:`\frac{1}{4}\left(1-\xi\right)\left(1-\eta\right)\left(1-\zeta^{2}\right)`
     - (:math:`-\frac{1}{4}\left(\eta-1\right)\left(\zeta^{2}-1\right)`, :math:`-\frac{1}{4}\left(\xi-1\right)\left(\zeta^{2}-1\right)`, :math:`3`)
   * - 14
     - (:math:`1`, :math:`-1`, :math:`0`)
     - :math:`\frac{1}{4}\left(1+\xi\right)\left(1-\eta\right)\left(1-\zeta^{2}\right)`
     - (:math:`\frac{1}{4}\left(\eta-1\right)\left(\zeta^{2}-1\right)`, :math:`\frac{1}{4}\left(\xi+1\right)\left(\zeta^{2}-1\right)`, :math:`3`)
   * - 15
     - (:math:`1`, :math:`1`, :math:`0`)
     - :math:`\frac{1}{4}\left(1+\xi\right)\left(1+\eta\right)\left(1-\zeta^{2}\right)`
     - (:math:`-\frac{1}{4}\left(\eta+1\right)\left(\zeta^{2}-1\right)`, :math:`-\frac{1}{4}\left(\xi+1\right)\left(\zeta^{2}-1\right)`, :math:`3`)
   * - 16
     - (:math:`-1`, :math:`1`, :math:`0`)
     - :math:`\frac{1}{4}\left(1-\xi\right)\left(1+\eta\right)\left(1-\zeta^{2}\right)`
     - (:math:`\frac{1}{4}\left(\eta+1\right)\left(\zeta^{2}-1\right)`, :math:`\frac{1}{4}\left(\xi-1\right)\left(\zeta^{2}-1\right)`, :math:`3`)
   * - 17
     - (:math:`0`, :math:`-1`, :math:`1`)
     - :math:`\frac{1}{4}\left(1-\xi^{2}\right)\left(1-\eta\right)\left(1+\zeta\right)`
     - (:math:`\frac{1}{2}\xi\left(\eta-1\right)\left(\zeta+1\right)`, :math:`\frac{1}{4}\left(\xi^{2}-1\right)\left(\zeta+1\right)`, :math:`3`)
   * - 18
     - (:math:`1`, :math:`0`, :math:`1`)
     - :math:`\frac{1}{4}\left(1+\xi\right)\left(1-\eta^{2}\right)\left(1+\zeta\right)`
     - (:math:`-\frac{1}{4}\left(\eta^{2}-1\right)\left(\zeta+1\right)`, :math:`-\frac{1}{2}\eta\left(\xi+1\right)\left(\zeta+1\right)`, :math:`3`)
   * - 19
     - (:math:`0`, :math:`1`, :math:`1`)
     - :math:`\frac{1}{4}\left(1-\xi^{2}\right)\left(1+\eta\right)\left(1+\zeta\right)`
     - (:math:`-\frac{1}{2}\xi\left(\eta+1\right)\left(\zeta+1\right)`, :math:`-\frac{1}{4}\left(\xi^{2}-1\right)\left(\zeta+1\right)`, :math:`3`)
   * - 20
     - (:math:`-1`, :math:`0`, :math:`1`)
     - :math:`\frac{1}{4}\left(1-\xi\right)\left(1-\eta^{2}\right)\left(1+\zeta\right)`
     - (:math:`\frac{1}{4}\left(\eta^{2}-1\right)\left(\zeta+1\right)`, :math:`\frac{1}{2}\eta\left(\xi-1\right)\left(\zeta+1\right)`, :math:`3`)

.. list-table:: Gaussian quadrature points
   :align: center

   * - Coord. (:math:`\xi`, :math:`\eta`, :math:`\zeta`)
     - Weight
   * - (:math:`-\sqrt{\tfrac{3}{5}}`, :math:`-\sqrt{\tfrac{3}{5}}`, :math:`-\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{125}{729}`
   * - (:math:`-\sqrt{\tfrac{3}{5}}`, :math:`-\sqrt{\tfrac{3}{5}}`, :math:`0`)
     - :math:`\frac{200}{729}`
   * - (:math:`-\sqrt{\tfrac{3}{5}}`, :math:`-\sqrt{\tfrac{3}{5}}`, :math:`\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{125}{729}`
   * - (:math:`-\sqrt{\tfrac{3}{5}}`, :math:`0`, :math:`-\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{200}{729}`
   * - (:math:`-\sqrt{\tfrac{3}{5}}`, :math:`0`, :math:`0`)
     - :math:`\frac{320}{729}`
   * - (:math:`-\sqrt{\tfrac{3}{5}}`, :math:`0`, :math:`\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{200}{729}`
   * - (:math:`-\sqrt{\tfrac{3}{5}}`, :math:`\sqrt{\tfrac{3}{5}}`, :math:`-\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{125}{729}`
   * - (:math:`-\sqrt{\tfrac{3}{5}}`, :math:`\sqrt{\tfrac{3}{5}}`, :math:`0`)
     - :math:`\frac{200}{729}`
   * - (:math:`-\sqrt{\tfrac{3}{5}}`, :math:`\sqrt{\tfrac{3}{5}}`, :math:`\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{125}{729}`
   * - (:math:`0`, :math:`-\sqrt{\tfrac{3}{5}}`, :math:`-\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{200}{729}`
   * - (:math:`0`, :math:`-\sqrt{\tfrac{3}{5}}`, :math:`0`)
     - :math:`\frac{320}{729}`
   * - (:math:`0`, :math:`-\sqrt{\tfrac{3}{5}}`, :math:`\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{200}{729}`
   * - (:math:`0`, :math:`0`, :math:`-\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{320}{729}`
   * - (:math:`0`, :math:`0`, :math:`0`)
     - :math:`\frac{512}{729}`
   * - (:math:`0`, :math:`0`, :math:`\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{320}{729}`
   * - (:math:`0`, :math:`\sqrt{\tfrac{3}{5}}`, :math:`-\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{200}{729}`
   * - (:math:`0`, :math:`\sqrt{\tfrac{3}{5}}`, :math:`0`)
     - :math:`\frac{320}{729}`
   * - (:math:`0`, :math:`\sqrt{\tfrac{3}{5}}`, :math:`\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{200}{729}`
   * - (:math:`\sqrt{\tfrac{3}{5}}`, :math:`-\sqrt{\tfrac{3}{5}}`, :math:`-\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{125}{729}`
   * - (:math:`\sqrt{\tfrac{3}{5}}`, :math:`-\sqrt{\tfrac{3}{5}}`, :math:`0`)
     - :math:`\frac{200}{729}`
   * - (:math:`\sqrt{\tfrac{3}{5}}`, :math:`-\sqrt{\tfrac{3}{5}}`, :math:`\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{125}{729}`
   * - (:math:`\sqrt{\tfrac{3}{5}}`, :math:`0`, :math:`-\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{200}{729}`
   * - (:math:`\sqrt{\tfrac{3}{5}}`, :math:`0`, :math:`0`)
     - :math:`\frac{320}{729}`
   * - (:math:`\sqrt{\tfrac{3}{5}}`, :math:`0`, :math:`\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{200}{729}`
   * - (:math:`\sqrt{\tfrac{3}{5}}`, :math:`\sqrt{\tfrac{3}{5}}`, :math:`-\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{125}{729}`
   * - (:math:`\sqrt{\tfrac{3}{5}}`, :math:`\sqrt{\tfrac{3}{5}}`, :math:`0`)
     - :math:`\frac{200}{729}`
   * - (:math:`\sqrt{\tfrac{3}{5}}`, :math:`\sqrt{\tfrac{3}{5}}`, :math:`\sqrt{\tfrac{3}{5}}`)
     - :math:`\frac{125}{729}`

Pentahedron 15
''''''''''''''

.. list-table:: Elements properties
   :header-rows: 1

   * - Node (:math:`i`)
     - Coord. (:math:`\xi`, :math:`\eta`, :math:`\zeta`)
     - Shape function (:math:`N_i`)
     - Derivative (:math:`\frac{\partial N_i}{\partial \xi}`, :math:`\frac{\partial N_i}{\partial \eta}`, :math:`\frac{\partial N_i}{\partial \zeta}`)
   * - 1
     - (:math:`-1`, :math:`1`, :math:`0`)
     - :math:`\frac{1}{2}\eta\left(1-\xi\right)\left(2\eta-2-\xi\right)`
     - (:math:`\frac{1}{2}\eta\left(2\xi-2\eta+1\right)`, :math:`-\frac{1}{2}\left(\xi-1\right)\left(4\eta-\xi-2\right)`, :math:`3`)
   * - 2
     - (:math:`-1`, :math:`0`, :math:`1`)
     - :math:`\frac{1}{2}\zeta\left(1-\xi\right)\left(2\zeta-2-\xi\right)`
     - (:math:`\frac{1}{2}\zeta\left(2\xi-2\zeta+1\right)`, :math:`0.0`, :math:`3`)
   * - 3
     - (:math:`-1`, :math:`0`, :math:`0`)
     - :math:`\frac{1}{2}\left(\xi-1\right)\left(1-\eta-\zeta\right)\left(\xi+2\eta+2\zeta\right)`
     - (:math:`-\frac{1}{2}\left(2\xi+2\eta+2\zeta-1\right)\left(\eta+\zeta-1\right)`, :math:`-\frac{1}{2}\left(\xi-1\right)\left(4\eta+\xi+2\left(2\zeta-1\right)\right)`, :math:`3`)
   * - 4
     - (:math:`1`, :math:`1`, :math:`0`)
     - :math:`\frac{1}{2}\eta\left(1+\xi\right)\left(2\eta-2+\xi\right)`
     - (:math:`\frac{1}{2}\eta\left(2\xi+2\eta-1\right)`, :math:`\frac{1}{2}\left(\xi+1\right)\left(4\eta+\xi-2\right)`, :math:`3`)
   * - 5
     - (:math:`1`, :math:`0`, :math:`1`)
     - :math:`\frac{1}{2}\zeta\left(1+\xi\right)\left(2\zeta-2+\xi\right)`
     - (:math:`\frac{1}{2}\zeta\left(2\xi+2\zeta-1\right)`, :math:`0.0`, :math:`3`)
   * - 6
     - (:math:`1`, :math:`0`, :math:`0`)
     - :math:`\frac{1}{2}\left(-\xi-1\right)\left(1-\eta-\zeta\right)\left(-\xi+2\eta+2\zeta\right)`
     - (:math:`-\frac{1}{2}\left(\eta+\zeta-1\right)\left(2\xi-2\eta-2\zeta+1\right)`, :math:`\frac{1}{2}\left(\xi+1\right)\left(4\eta-\xi+2\left(2\zeta-1\right)\right)`, :math:`3`)
   * - 7
     - (:math:`-1`, :math:`0.5`, :math:`0.5`)
     - :math:`2\eta\zeta\left(1-\xi\right)`
     - (:math:`-2\eta\zeta`, :math:`-2\left(\xi-1\right)\zeta`, :math:`3`)
   * - 8
     - (:math:`-1`, :math:`0`, :math:`0.5`)
     - :math:`2\zeta\left(1-\eta-\zeta\right)\left(1-\xi\right)`
     - (:math:`2\zeta\left(\eta+\zeta-1\right)`, :math:`2\zeta-\left(\xi-1\right)`, :math:`3`)
   * - 9
     - (:math:`-1`, :math:`0.5`, :math:`0`)
     - :math:`2\eta\left(1-\xi\right)\left(1-\eta-\zeta\right)`
     - (:math:`2\eta\left(\eta+\zeta-1\right)`, :math:`2\left(2\eta+\zeta-1\right)\left(\xi-1\right)`, :math:`3`)
   * - 10
     - (:math:`0`, :math:`1`, :math:`0`)
     - :math:`\eta\left(1-\xi^{2}\right)`
     - (:math:`-2\xi\eta`, :math:`-\left(\xi^{2}-1\right)`, :math:`3`)
   * - 11
     - (:math:`0`, :math:`0`, :math:`1`)
     - :math:`\zeta\left(1-\xi^{2}\right)`
     - (:math:`-2\xi\zeta`, :math:`0.0`, :math:`3`)
   * - 12
     - (:math:`0`, :math:`0`, :math:`0`)
     - :math:`\left(1-\xi^{2}\right)\left(1-\eta-\zeta\right)`
     - (:math:`2\xi\left(\eta+\zeta-1\right)`, :math:`\left(\xi^{2}-1\right)`, :math:`3`)
   * - 13
     - (:math:`1`, :math:`0.5`, :math:`0.5`)
     - :math:`2\eta\zeta\left(1+\xi\right)`
     - (:math:`2\eta\zeta`, :math:`2\zeta\left(\xi+1\right)`, :math:`3`)
   * - 14
     - (:math:`1`, :math:`0`, :math:`0.5`)
     - :math:`2\zeta\left(1+\xi\right)\left(1-\eta-\zeta\right)`
     - (:math:`-2\zeta\left(\eta+\zeta-1\right)`, :math:`-2\zeta\left(\xi+1\right)`, :math:`3`)
   * - 15
     - (:math:`1`, :math:`0.5`, :math:`0`)
     - :math:`2\eta\left(1+\xi\right)\left(1-\eta-\zeta\right)`
     - (:math:`-2\eta\left(\eta+\zeta-1\right)`, :math:`-2\left(2\eta+\zeta-1\right)\left(\xi+1\right)`, :math:`3`)

.. list-table:: Gaussian quadrature points
   :align: center

   * - Coord. (:math:`\xi`, :math:`\eta`, :math:`\zeta`)
     - Weight
   * - (:math:`-{\tfrac{1}{\sqrt{3}}}`, :math:`\tfrac{1}{3}`, :math:`\tfrac{1}{3}`)
     - -:math:`\frac{27}{96}`
   * - (:math:`-{\tfrac{1}{\sqrt{3}}}`, :math:`0.6`, :math:`0.2`)
     - :math:`\frac{25}{96}`
   * - (:math:`-{\tfrac{1}{\sqrt{3}}}`, :math:`0.2`, :math:`0.6`)
     - :math:`\frac{25}{96}`
   * - (:math:`-{\tfrac{1}{\sqrt{3}}}`, :math:`0.2`, :math:`0.2`)
     - :math:`\frac{25}{96}`
   * - (:math:`{\tfrac{1}{\sqrt{3}}}`, :math:`\tfrac{1}{3}`, :math:`\tfrac{1}{3}`)
     - -:math:`\frac{27}{96}`
   * - (:math:`{\tfrac{1}{\sqrt{3}}}`, :math:`0.6`, :math:`0.2`)
     - :math:`\frac{25}{96}`
   * - (:math:`{\tfrac{1}{\sqrt{3}}}`, :math:`0.2`, :math:`0.6`)
     - :math:`\frac{25}{96}`
   * - (:math:`{\tfrac{1}{\sqrt{3}}}`, :math:`0.2`, :math:`0.2`)
     - :math:`\frac{25}{96}`
