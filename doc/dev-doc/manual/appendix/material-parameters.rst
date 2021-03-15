.. _app-material-parameters:

Material Parameters
===================

Linear elastic isotropic
------------------------

Keyword: :ref:`elastic <sect-smm-linear-elastic-isotropic>`

Parameters:

- ``rho``: (*Real*) Density
- ``E``: (*Real*) Young's modulus
- ``nu``: (*Real*) Poisson's ratio
- ``Plane_stress``: (*bool*) Plane stress simplification (only 2D problems)


Linear elastic anisotropic
--------------------------

Keyword: :ref:`elastic_anisotropic <sect-smm-linear-elastic-anisotropic>`

Parameters:

- ``rho``: (*Real*) Density
- ``n1``: (*Vector<Real>*) Direction of main material axis
- ``n2``: (*Vector<Real>*) Direction of second material axis
- ``n3``: (*Vector<Real>*) Direction of third material axis
- ``C..``: (*Real*) Coefficient ij of material tensor C (all the 36 values in
  Voigt notation can be entered)
- ``alpha``: (*Real*) Viscous propertion (default is 0)


Linear elastic anisotropic
------------------------

Keyword: :ref:`elastic_orthotropic <sect-smm-linear-elastic-orthotropic>`

Parameters:

- ``rho``: (*Real*) Density
- ``n1``: (*Vector<Real>*) Direction of main material axis
- ``n2``: (*Vector<Real>*) Direction of second material axis (if applicable)
- ``n3``: (*Vector<Real>*) Direction of third material axis (if applicable)
- ``E1``: (*Real*) Young's modulus (n1)
- ``E2``: (*Real*) Young's modulus (n2)
- ``E3``: (*Real*) Young's modulus (n3)
- ``nu1``: (*Real*) Poisson's ratio (n1)
- ``nu2``: (*Real*) Poisson's ratio (n2)
- ``nu3``: (*Real*) Poisson's ratio (n3)
- ``G12``: (*Real*) Shear modulus (12)
- ``G13``: (*Real*) Shear modulus (13)
- ``G23``: (*Real*) Shear modulus (23)


Neohookean (finite strains)
---------------------------

Keyword: :ref:`neohookean <sect-smm-cl-neohookean>`

Parameters:

- ``rho``: (*Real*) Density
- ``E``: (*Real*) Young's modulus
- ``nu``: (*Real*) Poisson's ratio
- ``Plane_stress``: (*bool*) Plane stress simplification (only 2D problems)


Standard linear solid
---------------------

Keyword: :ref:`sls_deviatoric <sect-smm-cl-sls>`

Parameters:

- ``rho``: (*Real*) Density
- ``E``: (*Real*) Young's modulus
- ``nu``: (*Real*) Poisson's ratio
- ``Plane_stress``: (*bool*) Plane stress simplification (only 2D problems)
- ``Eta``: (*Real*) Viscosity
- ``Ev``: (*Real*) Stiffness of viscous element


Elasto-plastic linear isotropic hardening
-----------------------------------------

Keyword: :ref:`plastic_linear_isotropic_hardening <sect-smm-cl-plastic>`

Parameters:

- ``rho``: (*Real*) Density
- ``E``: (*Real*) Young's modulus
- ``nu``: (*Real*) Poisson's ratio
- ``h``: (*Real*) Hardening modulus
- ``sigma_y``: (*Real*) Yield stress


Marigo
------

Keyword: :ref:`marigo <sect-smm-cl-damage-marigo>`

Parameters:

- ``rho``: (*Real*) Density
- ``E``: (*Real*) Young's modulus
- ``nu``: (*Real*) Poisson's ratio
- ``Plane_stress``: (*bool*) Plane stress simplification (only 2D problems)
- ``Yd``: (*Random*) Hardening modulus
- ``Sd``: (*Real*) Damage energy


Mazars
------

Keyword: :ref:`mazars <sect-smm-cl-damage-mazars>`

Parameters:

- ``rho``: (*Real*) Density
- ``E``: (*Real*) Young's modulus
- ``nu``: (*Real*) Poisson's ratio
- ``At``: (*Real*) Traction post-peak asymptotic value
- ``Bt``: (*Real*) Traction decay shape
- ``Ac``: (*Real*) Compression post-peak asymptotic value
- ``Bc``: (*Real*) Compression decay shape
- ``K0``: (*Real*) Damage threshold
- ``beta``: (*Real*) Shear parameter
