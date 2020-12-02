.. _sect-smm-cl:

Constitutive Laws
-----------------

In order to compute an element’s response to deformation, one needs to
use an appropriate constitutive relationship. The constitutive law is
used to compute the element’s stresses from the element’s strains.

In the finite-element discretization, the constitutive formulation is
applied to every quadrature point of each element. When the implicit
formulation is used, the tangent matrix has to be computed.

| The chosen materials for the simulation have to be specified in the
  mesh file or, as an alternative, they can be assigned using the at
  ``element_material`` vector. For
  every material assigned to the problem one has to specify the material
  characteristics (constitutive behavior and material properties) using
  the text input file (see `[sect:io:material] <#sect:io:material>`__).
| In order to conveniently store values at each quadrature in a material point
  ``Akantu`` provides a special data structure, the at :cpp:class:`InternalField
  <akantu::InternalField>`. The internal fields are inheriting from the at
  :cpp:class:`ElementTypeMapArray <akantu::ElementTypeMapArray>`. Furthermore,
  it provides several functions for initialization, auto-resizing and auto
  removal of quadrature points.

Sometimes it is also desired to generate random distributions of
internal parameters. An example might be the critical stress at which
the material fails. To generate such a field, in the text input file, a
random quantity needs be added to the base value:

All parameters are real numbers. For the uniform distribution, minimum
and maximum values have to be specified. Random parameters are defined
as a :math:`base` value to which we add a random number that follows the
chosen distribution.

The
`Uniform <http://en.wikipedia.org/wiki/Uniform_distribution_(continuous)>`__
distribution is gives a random values between in :math:`[min, max)`. The
`Weibull <http://en.wikipedia.org/wiki/Weibull_distribution>`__
distribution is characterized by the following cumulative distribution
function:

.. math:: F(x) = 1- e^{-\left({x/\lambda}\right)^m}

which depends on :math:`m` and :math:`\lambda`, which are the shape
parameter and the scale parameter. These random distributions are
different each time the code is executed. In order to obtain always the
same one, it possible to manually set the *seed* that is the number from
which these pseudo-random distributions are created. This can be done by
adding the following line to the input file *outside* the material
parameters environments:

.. code-block::

   seed = 1.0

where the value 1.0 can be substituted with any number. Currently
``Akantu`` can reproduce always the same distribution when the seed is
specified *only* in serial. The value of the *seed* can be also
specified directly in the code (for instance in the main file) with the
command:

.. code-block::

   RandGenerator<Real>::seed(1.0)

The same command, with empty brackets, can be used to check the value of
the *seed* used in the simulation.

The following sections describe the constitutive models implemented in
``Akantu``. In Appendix `7 <#app:material-parameters>`__ a summary of
the parameters for all materials of ``Akantu`` is provided.

Elasticity
``````````

The elastic law is a commonly used constitutive relationship that can be
used for a wide range of engineering materials (*e.g.*, metals,
concrete, rock, wood, glass, rubber, etc.) provided that the strains
remain small (*i.e.*, small deformation and stress lower than yield
strength).

The elastic laws are often expressed as
:math:`\boldsymbol{\sigma} =
\boldsymbol{C}:\boldsymbol{\varepsilon}` with
where :math:`\boldsymbol{\sigma}` is the Cauchy stress
tensor, :math:`\boldsymbol{\varepsilon}` represents the
infinitesimal strain tensor and :math:`\boldsymbol{C}` is
the elastic modulus tensor.

.. _sect-smm-linear-elastic-isotropic:

Linear isotropic
''''''''''''''''

The linear isotropic elastic behavior is described by Hooke’s law, which
states that the stress is linearly proportional to the applied strain
(material behaves like an ideal spring), as illustrated in
Figure `[fig:smm:cl:elastic] <#fig:smm:cl:elastic>`__.

The equation that relates the strains to the displacements is: point)
from the displacements as follows:

.. math::

   \label{eqn:smm:strain_inf}
     \boldsymbol{\varepsilon} =
     \frac{1}{2} \left[ \nabla_0 \boldsymbol{u}+\nabla_0 \boldsymbol{u}^T \right]

where :math:`\boldsymbol{\varepsilon}` represents the
infinitesimal strain tensor,
:math:`\nabla_{0}\boldsymbol{u}` the displacement gradient
tensor according to the initial configuration. The constitutive equation
for isotropic homogeneous media can be expressed as:

.. math::

   \label{eqn:smm:material:constitutive_elastic}
     \boldsymbol{\sigma } =\lambda\mathrm{tr}(\boldsymbol{\varepsilon})\boldsymbol{I}+2 \mu\boldsymbol{\varepsilon}

where :math:`\boldsymbol{\sigma}` is the Cauchy stress
tensor (:math:`\lambda` and :math:`\mu` are the the first and second
Lame’s coefficients).

In Voigt notation this correspond to

.. math::

   \begin{aligned}
     \left[\begin{array}{c}
         \sigma_{11}\\
         \sigma_{22}\\
         \sigma_{33}\\
         \sigma_{23}\\
         \sigma_{13}\\
         \sigma_{12}\\
       \end{array}\right]
     &= \frac{E}{(1+\nu)(1-2\nu)}\left[
       \begin{array}{cccccc}
         1-\nu & \nu   & \nu   & 0 & 0 & 0\\
         \nu   & 1-\nu & \nu   & 0 & 0 & 0\\
         \nu   & \nu   & 1-\nu & 0 & 0 & 0\\
         0     &  0    &  0    & \frac{1-2\nu}{2} & 0 & 0 \\
         0     &  0    &  0    & 0 & \frac{1-2\nu}{2} & 0 \\
         0     &  0    &  0    & 0 & 0 & \frac{1-2\nu}{2} \\
       \end{array}\right]
     \left[\begin{array}{c}
         \varepsilon_{11}\\
         \varepsilon_{22}\\
         \varepsilon_{33}\\
         2\varepsilon_{23}\\
         2\varepsilon_{13}\\
         2\varepsilon_{12}\\
       \end{array}\right]\end{aligned}

.. _sect-smm-linear-elastic-anisotropic:

Linear anisotropic
''''''''''''''''''

This formulation is not sufficient to represent all elastic material
behavior. Some materials have characteristic orientation that have to be
taken into account. To represent this anisotropy a more general
stress-strain law has to be used. For this we define the elastic modulus
tensor as follow:

.. math::

   \begin{aligned}
     \left[\begin{array}{c}
         \sigma_{11}\\
         \sigma_{22}\\
         \sigma_{33}\\
         \sigma_{23}\\
         \sigma_{13}\\
         \sigma_{12}\\
       \end{array}\right]
     &= \left[
       \begin{array}{cccccc}
         c_{11} & c_{12} & c_{13} & c_{14} & c_{15} & c_{16}\\
         c_{21} & c_{22} & c_{23} & c_{24} & c_{25} & c_{26}\\
         c_{31} & c_{32} & c_{33} & c_{34} & c_{35} & c_{36}\\
         c_{41} & c_{42} & c_{43} & c_{44} & c_{45} & c_{46}\\
         c_{51} & c_{52} & c_{53} & c_{54} & c_{55} & c_{56}\\
         c_{61} & c_{62} & c_{63} & c_{64} & c_{65} & c_{66}\\
       \end{array}\right]
     \left[\begin{array}{c}
         \varepsilon_{11}\\
         \varepsilon_{22}\\
         \varepsilon_{33}\\
         2\varepsilon_{23}\\
         2\varepsilon_{13}\\
         2\varepsilon_{12}\\
       \end{array}\right]\end{aligned}

To simplify the writing of input files the :math:`\boldsymbol{C}` tensor
is expressed in the material basis. And this basis as to be given too.
This basis :math:`\Omega_{{\mathrm{mat}}}
= \{\boldsymbol{n_1}, \boldsymbol{n_2}, \boldsymbol{n_3}\}`
is used to define the rotation :math:`R_{ij} =
\boldsymbol{n_j} . \boldsymbol{e_i}`. And
:math:`\boldsymbol{C}` can be rotated in the global basis
:math:`\Omega
= \{\boldsymbol{e_1}, \boldsymbol{e_2}, \boldsymbol{e_3}\}`
as follow:

.. math::

   \begin{aligned}
   \boldsymbol{C}_{\Omega} &= \boldsymbol{R}_1 \boldsymbol{C}_{\Omega_{{\mathrm{mat}}}} \boldsymbol{R}_2\\
   \boldsymbol{R}_1  &= \left[
     \begin{array}{cccccc}
       R_{11} R_{11} & R_{12} R_{12} & R_{13} R_{13} & R_{12} R_{13} & R_{11} R_{13} & R_{11} R_{12}\\
       R_{21} R_{21} & R_{22} R_{22} & R_{23} R_{23} & R_{22} R_{23} & R_{21} R_{23} & R_{21} R_{22}\\
       R_{31} R_{31} & R_{32} R_{32} & R_{33} R_{33} & R_{32} R_{33} & R_{31} R_{33} & R_{31} R_{32}\\
       R_{21} R_{31} & R_{22} R_{32} & R_{23} R_{33} & R_{22} R_{33} & R_{21} R_{33} & R_{21} R_{32}\\
       R_{11} R_{31} & R_{12} R_{32} & R_{13} R_{33} & R_{12} R_{33} & R_{11} R_{33} & R_{11} R_{32}\\
       R_{11} R_{21} & R_{12} R_{22} & R_{13} R_{23} & R_{12} R_{23} & R_{11} R_{23} & R_{11} R_{22}\\
     \end{array}\right]\\
   \boldsymbol{R}_2  &= \left[
     \begin{array}{cccccc}
       R_{11} R_{11} & R_{21} R_{21} & R_{31} R_{31} & R_{21} R_{31} & R_{11} R_{31} & R_{11} R_{21}\\
       R_{12} R_{12} & R_{22} R_{22} & R_{32} R_{32} & R_{22} R_{32} & R_{12} R_{32} & R_{12} R_{22}\\
       R_{13} R_{13} & R_{23} R_{23} & R_{33} R_{33} & R_{23} R_{33} & R_{13} R_{33} & R_{13} R_{23}\\
       R_{12} R_{13} & R_{22} R_{23} & R_{32} R_{33} & R_{22} R_{33} & R_{12} R_{33} & R_{12} R_{23}\\
       R_{11} R_{13} & R_{21} R_{23} & R_{31} R_{33} & R_{21} R_{33} & R_{11} R_{33} & R_{11} R_{23}\\
       R_{11} R_{12} & R_{21} R_{22} & R_{31} R_{32} & R_{21} R_{32} & R_{11} R_{32} & R_{11} R_{22}\\
     \end{array}\right]\\\end{aligned}

.. _sect-smm-linear-elastic-orthotropic:

Linear orthotropic
''''''''''''''''''

A particular case of anisotropy is when the material basis is orthogonal
in which case the elastic modulus tensor can be simplified and rewritten
in terms of 9 independents material parameters.

.. math::

   \begin{aligned}
     \left[\begin{array}{c}
         \sigma_{11}\\
         \sigma_{22}\\
         \sigma_{33}\\
         \sigma_{23}\\
         \sigma_{13}\\
         \sigma_{12}\\
       \end{array}\right]
     &= \left[
       \begin{array}{cccccc}
         c_{11} & c_{12} & c_{13} &   0   &   0   &   0  \\
               & c_{22} & c_{23} &   0   &   0   &   0  \\
               &       & c_{33} &   0   &   0   &   0  \\
               &       &       & c_{44} &   0   &   0  \\
               &  \multicolumn{2}{l}{\text{sym.}}       &       & c_{55} &   0  \\
               &       &       &       &       & c_{66}\\
       \end{array}\right]
     \left[\begin{array}{c}
         \varepsilon_{11}\\
         \varepsilon_{22}\\
         \varepsilon_{33}\\
         2\varepsilon_{23}\\
         2\varepsilon_{13}\\
         2\varepsilon_{12}\\
       \end{array}\right]\end{aligned}

.. math::

   \begin{aligned}
     c_{11} &= E_1 (1 - \nu_{23}\nu_{32})\Gamma \qquad c_{22} = E_2 (1 - \nu_{13}\nu_{31})\Gamma \qquad c_{33} = E_3 (1 - \nu_{12}\nu_{21})\Gamma\\
     c_{12} &= E_1 (\nu_{21} - \nu_{31}\nu_{23})\Gamma = E_2 (\nu_{12} - \nu_{32}\nu_{13})\Gamma\\
     c_{13} &= E_1 (\nu_{31} - \nu_{21}\nu_{32})\Gamma = E_2 (\nu_{13} - \nu_{21}\nu_{23})\Gamma\\
     c_{23} &= E_2 (\nu_{32} - \nu_{12}\nu_{31})\Gamma = E_3 (\nu_{23} - \nu_{21}\nu_{13})\Gamma\\
     c_{44} &= \mu_{23} \qquad  c_{55} = \mu_{13} \qquad  c_{66} = \mu_{12} \\
     \Gamma &= \frac{1}{1 - \nu_{12} \nu_{21} - \nu_{13} \nu_{31} - \nu_{32} \nu_{23} - 2 \nu_{21} \nu_{32} \nu_{13}}\end{aligned}

The Poisson ratios follow the rule
:math:`\nu_{ij} = \nu_{ji} E_i / E_j`.

.. _sect-smm-cl-neohookean:

Neo-Hookean
'''''''''''

The hyperelastic Neo-Hookean constitutive law results from an extension
of the linear elastic relationship (Hooke’s Law) for large deformation.
Thus, the model predicts nonlinear stress-strain behavior for bodies
undergoing large deformations.

.. figure:: figures/stress_strain_neo.pdf
   :alt: Neo-hookean Stress-strain curve.
   :name: fig:smm:cl:neo_hookean
   :width: 40.0%

   Neo-hookean Stress-strain curve.

As illustrated in Figure `4.6 <#fig:smm:cl:neo_hookean>`__, the behavior
is initially linear and the mechanical behavior is very close to the
corresponding linear elastic material. This constitutive relationship,
which accounts for compressibility, is a modified version of the one
proposed by Ronald Rivlin :cite:`Belytschko:2000`.

The strain energy stored in the material is given by:

.. math::

   \label{eqn:smm:constitutive:neohookean_potential}
     \Psi(\boldsymbol{C}) = \frac{1}{2}\lambda_0\left(\ln J\right)^2-\mu_0\ln J+\frac{1}{2}
     \mu_0\left(\mathrm{tr}(\boldsymbol{C})-3\right)

where :math:`\lambda_0` and :math:`\mu_0` are, respectively, Lamé’s
first parameter and the shear modulus at the initial configuration.
:math:`J` is the jacobian of the deformation gradient
(:math:`\boldsymbol{F}=\nabla_{\!\!\boldsymbol{X}}\boldsymbol{x}`):
:math:`J=\text{det}(\boldsymbol{F})`. Finally
:math:`\boldsymbol{C}` is the right Cauchy-Green
deformation tensor.

Since this kind of material is used for large deformation problems, a
finite deformation framework should be used. Therefore, the Cauchy
stress (:math:`\boldsymbol{\sigma}`) should be computed
through the second Piola-Kirchhoff stress tensor
:math:`\boldsymbol{S}`:

.. math:: \boldsymbol{\sigma } = \frac{1}{J}\boldsymbol{F}\boldsymbol{S}\boldsymbol{F}^T

Finally the second Piola-Kirchhoff stress tensor is given by:

.. math::

   \boldsymbol{S}  = 2\frac{\partial\Psi}{\partial\boldsymbol{C}} = \lambda_0\ln J
     \boldsymbol{C}^{-1}+\mu_0\left(\boldsymbol{I}-\boldsymbol{C}^{-1}\right)

The parameters to indicate in the material file are the same as those
for the elastic case: ``E`` (Young’s modulus), ``nu`` (Poisson’s ratio).

.. _sect-smm-cl-sls:

Visco-Elasticity
''''''''''''''''

Visco-elasticity is characterized by strain rate dependent behavior.
Moreover, when such a material undergoes a deformation it dissipates
energy. This dissipation results in a hysteresis loop in the
stress-strain curve at every loading cycle (see
Figure `[fig:smm:cl:visco-elastic:hyst] <#fig:smm:cl:visco-elastic:hyst>`__).
In principle, it can be applied to many materials, since all materials
exhibit a visco-elastic behavior if subjected to particular conditions
(such as high temperatures).

The standard rheological linear solid model (see Sections 10.2 and 10.3
of :cite:`simo92`) has been implemented in ``Akantu``. This
model results from the combination of a spring mounted in parallel with
a spring and a dashpot connected in series, as illustrated in
Figure `[fig:smm:cl:visco-elastic:model] <#fig:smm:cl:visco-elastic:model>`__.
The advantage of this model is that it allows to account for creep or
stress relaxation. The equation that relates the stress to the strain is
(in 1D):

.. math:: \frac{d\varepsilon(t)}{dt} = \left ( E + E_V \right ) ^ {-1} \cdot \left [ \frac{d\sigma(t)}{dt} + \frac{E_V}{\eta}\sigma(t) - \frac{EE_V}{\eta}\varepsilon(t) \right ]

where :math:`\eta` is the viscosity. The equilibrium condition is unique and is
attained in the limit, as :math:`t \to \infty`. At this stage, the response is
elastic and depends on the Young’s modulus :math:`E`. The mandatory parameters
for the material file are the following: ``rho`` (density), ``E`` (Young’s
modulus), ``nu`` (Poisson’s ratio), ``Plane_Stress`` (if set to zero plane
strain, otherwise plane stress), ``eta`` (dashpot viscosity) and ``Ev``
(stiffness of the viscous element).

Note that the current standard linear solid model is applied only on the
deviatoric part of the strain tensor. The spheric part of the strain
tensor affects the stress tensor like an linear elastic material.

.. _sect-smm-cl-plastic:

Small-Deformation Plasticity
''''''''''''''''''''''''''''

The small-deformation plasticity is a simple plasticity material
formulation which accounts for the additive decomposition of strain into
elastic and plastic strain components. This formulation is applicable to
infinitesimal deformation where the additive decomposition of the strain
is a valid approximation. In this formulation, plastic strain is a
shearing process where hydrostatic stress has no contribution to
plasticity and consequently plasticity does not lead to volume change.
Figure `4.7 <#fig:smm:cl:Lin-strain-hard>`__ shows the linear strain
hardening elasto-plastic behavior according to the additive
decomposition of strain into the elastic and plastic parts in
infinitesimal deformation as

.. math::

   \boldsymbol{\varepsilon} &= \boldsymbol{\varepsilon}^e +\boldsymbol{\varepsilon}^p\\
   \boldsymbol{\sigma} &= 2G(\boldsymbol{\varepsilon}^e) + \lambda  \mathrm{tr}(\boldsymbol{\varepsilon}^e)\boldsymbol{I}


In this class, the von Mises yield criterion is used. In the von Mises
yield criterion, the yield is independent of the hydrostatic stress.
Other yielding criteria such as Tresca and Gurson can be easily
implemented in this class as well.

In the von Mises yield criterion, the hydrostatic stresses have no
effect on the plasticity and consequently the yielding occurs when a
critical elastic shear energy is achieved.

.. math::

   \label{eqn:smm:constitutive:von Mises}
     f = \sigma_{{\mathrm{eff}}} - \sigma_y = \left(\frac{3}{2} {\boldsymbol{\sigma}}^{{\mathrm{tr}}} : {\boldsymbol{\sigma}}^{{\mathrm{tr}}}\right)^\frac{1}{2}-\sigma_y (\boldsymbol{\varepsilon}^p)

.. math::

   \label{eqn:smm:constitutive:yielding}
     f < 0 \quad \textrm{Elastic deformation,} \qquad f = 0 \quad  \textrm{Plastic deformation}

where :math:`\sigma_y` is the yield strength of the material which can
be function of plastic strain in case of hardening type of materials and
:math:`{\boldsymbol{\sigma}}^{{\mathrm{tr}}}` is the
deviatoric part of stress given by

.. math::

   \label{eqn:smm:constitutive:deviatoric stress}
     {\boldsymbol{\sigma}}^{{\mathrm{tr}}}=\boldsymbol{\sigma} - \frac{1}{3} \mathrm{tr}(\boldsymbol{\sigma}) \boldsymbol{I}

After yielding :math:`(f = 0)`, the normality hypothesis of plasticity
determines the direction of plastic flow which is normal to the tangent
to the yielding surface at the load point. Then, the tensorial form of
the plastic constitutive equation using the von Mises yielding criterion
(see equation 4.34) may be written as

.. math::

   \label{eqn:smm:constitutive:plastic contitutive equation}
     \Delta {\boldsymbol{\varepsilon}}^p = \Delta p \frac {\partial{f}}{\partial{\boldsymbol{\sigma}}}=\frac{3}{2} \Delta p \frac{{\boldsymbol{\sigma}}^{{\mathrm{tr}}}}{\sigma_{{\mathrm{eff}}}}

In these expressions, the direction of the plastic strain increment (or
equivalently, plastic strain rate) is given by
:math:`\frac{{\boldsymbol{\sigma}}^{{\mathrm{tr}}}}{\sigma_{{\mathrm{eff}}}}`
while the magnitude is defined by the plastic multiplier
:math:`\Delta p`. This can be obtained using the *consistency condition*
which impose the requirement for the load point to remain on the
yielding surface in the plastic regime.

Here, we summarize the implementation procedures for the
small-deformation plasticity with linear isotropic hardening:

#. Compute the trial stress:

   .. math:: {\boldsymbol{\sigma}}^{{\mathrm{tr}}} = {\boldsymbol{\sigma}}_t + 2G\Delta \boldsymbol{\varepsilon} + \lambda \mathrm{tr}(\Delta \boldsymbol{\varepsilon})\boldsymbol{I}

#. Check the Yielding criteria:

   .. math:: f = (\frac{3}{2} {\boldsymbol{\sigma}}^{{\mathrm{tr}}} : {\boldsymbol{\sigma}}^{{\mathrm{tr}}})^{1/2}-\sigma_y (\boldsymbol{\varepsilon}^p)

#. Compute the Plastic multiplier:

   .. math::

      \begin{aligned}
          d \Delta p &= \frac{\sigma^{tr}_{eff} - 3G \Delta P^{(k)}- \sigma_y^{(k)}}{3G + h}\\
          \Delta p^{(k+1)} &= \Delta p^{(k)}+ d\Delta p\\
          \sigma_y^{(k+1)} &= (\sigma_y)_t+ h\Delta p
        \end{aligned}

#. Compute the plastic strain increment:

   .. math:: \Delta {\boldsymbol{\varepsilon}}^p = \frac{3}{2} \Delta p \frac{{\boldsymbol{\sigma}}^{{\mathrm{tr}}}}{\sigma_{{\mathrm{eff}}}}

#. Compute the stress increment:

   .. math:: {\Delta \boldsymbol{\sigma}} = 2G(\Delta \boldsymbol{\varepsilon}-\Delta \boldsymbol{\varepsilon}^p) + \lambda  \mathrm{tr}(\Delta \boldsymbol{\varepsilon}-\Delta \boldsymbol{\varepsilon}^p)\boldsymbol{I}

#. Update the variables:

   .. math::

      \begin{aligned}
          {\boldsymbol{\varepsilon^p}} &= {\boldsymbol{\varepsilon}}^p_t+{\Delta {\boldsymbol{\varepsilon}}^p}\\
          {\boldsymbol{\sigma}} &= {\boldsymbol{\sigma}}_t+{\Delta \boldsymbol{\sigma}}
        \end{aligned}

We use an implicit integration technique called *the radial return method* to
obtain the plastic multiplier. This method has the advantage of being
unconditionally stable, however, the accuracy remains dependent on the step
size. The plastic parameters to indicate in the material file are:
:math:`\sigma_y` (Yield stress) and ``h`` (Hardening modulus). In addition, the
elastic parameters need to be defined as previously mentioned: ``E`` (Young’s
modulus), ``nu`` (Poisson’s ratio).

Damage
``````

In the simplified case of a linear elastic and brittle material,
isotropic damage can be represented by a scalar variable :math:`d`,
which varies from :math:`0` to :math:`1` for no damage to fully broken
material respectively. The stress-strain relationship then becomes:

.. math:: \boldsymbol{\sigma} = (1-d)\, \boldsymbol{C}:\boldsymbol{\varepsilon}

where :math:`\boldsymbol{\sigma}`,
:math:`\boldsymbol{\varepsilon}` are the Cauchy stress and
strain tensors, and :math:`\boldsymbol{C}` is the elastic
stiffness tensor. This formulation relies on the definition of an
evolution law for the damage variable. In ``Akantu``, many possibilities
exist and they are listed below.

.. _sect-smm-cl-damage:

Marigo
''''''

This damage evolution law is energy based as defined by Marigo
:cite:`marigo81a, lemaitre96a`. It is an isotropic damage law.

.. math::

   \begin{aligned}
     Y &= \frac{1}{2}\boldsymbol{\varepsilon}:\boldsymbol{C}:\boldsymbol{\varepsilon}\\
     F &= Y - Y_d - S d\\
     d &= \left\{
       \begin{array}{l l}
         \mathrm{min}\left(\frac{Y-Y_d}{S},\;1\right) & \mathrm{if}\; F > 0\\
         \mathrm{unchanged} & \mathrm{otherwise}
       \end{array}
     \right.\end{aligned}

In this formulation, :math:`Y` is the strain energy release rate,
:math:`Y_d` the rupture criterion and :math:`S` the damage energy. The
non-local version of this damage evolution law is constructed by
averaging the energy :math:`Y`.

.. _sect-smm-cl-damage-mazars:

Mazars
''''''

This law introduced by Mazars :cite:`mazars84a` is a
behavioral model to represent damage evolution in concrete. This model
does not rely on the computation of the tangent stiffness, the damage is
directly evaluated from the strain.

The governing variable in this damage law is the equivalent strain
:math:`\varepsilon_{{\mathrm{eq}}} =
\sqrt{<\boldsymbol{\varepsilon}>_+:<\boldsymbol{\varepsilon}>_+}`,
with :math:`<.>_+` the positive part of the tensor. This part is defined
in the principal coordinates (I, II, III) as
:math:`\varepsilon_{{\mathrm{eq}}} =
\sqrt{<\boldsymbol{\varepsilon_I}>_+^2 + <\boldsymbol{\varepsilon_{II}}>_+^2 + <\boldsymbol{\varepsilon_{III}}>_+^2}`.
The damage is defined as:

.. math::

   \begin{aligned}
     D &= \alpha_t^\beta D_t + (1-\alpha_t)^\beta D_c\\
     D_t &= 1 - \frac{\kappa_0 (1- A_t)}{\varepsilon_{{\mathrm{eq}}}} - A_t \exp^{-B_t(\varepsilon_{{\mathrm{eq}}}-\kappa_0)}\\
     D_c &= 1 - \frac{\kappa_0 (1- A_c)}{\varepsilon_{{\mathrm{eq}}}} - A_c
     \exp^{-B_c(\varepsilon_{{\mathrm{eq}}}-\kappa_0)}\\
     \alpha_t &= \frac{\sum_{i=1}^3<\varepsilon_i>_+\varepsilon_{{\mathrm{nd}}\;i}}{\varepsilon_{{\mathrm{eq}}}^2}\end{aligned}

With :math:`\kappa_0` the damage threshold, :math:`A_t` and :math:`B_t`
the damage parameter in traction, :math:`A_c` and :math:`B_c` the damage
parameter in compression, :math:`\beta` is the shear parameter.
:math:`\alpha_t` is the coupling parameter between traction and
compression, the :math:`\varepsilon_i` are the eigenstrain and the
:math:`\varepsilon_{{\mathrm{nd}}\;i}` are the eigenvalues of the strain
if the material were undamaged.

The coefficients :math:`A` and :math:`B` are the post-peak asymptotic
value and the decay shape parameters.
