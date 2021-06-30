Adding a New Constitutive Law
-----------------------------

There are several constitutive laws in ``Akantu`` as described in the previous
Section :ref:`sect-smm-cl`. It is also possible to use a user-defined material
for the simulation. These materials are referred to as local materials since
they are local to the example of the user and not part of the ``Akantu``
library. To define a new local material, two files (``material_XXX.hh`` and
``material_XXX.cc``) have to be provided where ``XXX`` is the name of the new
material. The header file ``material_XXX.hh`` defines the interface of your
custom material. Its implementation is provided in the ``material_XXX.cc``. The
new law must inherit from the :cpp:class:`Material <akantu::Material>` class or
any other existing material class. It is therefore necessary to include the
interface of the parent material in the header file of your local material and
indicate the inheritance in the declaration of the class::

   auto & solver = model.getNonLinearSolver();
   solver.set("max_iterations", 1);
   solver.set("threshold", 1e-4);
   solver.set("convergence_type", SolveConvergenceCriteria::_residual);

   model.solveStep();


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
  example, the Lam\'{e} constants of elastic materials can be considered as such
  parameters.

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
