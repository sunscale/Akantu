/**
 * @file   shape_igfem_inline_impl.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  ShapeIGFEM inline implementation
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_memory.hh"
//#include "mesh.hh"
#include "shape_igfem.hh"
/* -------------------------------------------------------------------------- */

#if defined(AKANTU_IGFEM)

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
ShapeLagrange<_ek_igfem>::ShapeLagrange(const Mesh & mesh,
				   const ID & id,
				   const MemoryID & memory_id) :
  ShapeFunctions(mesh, id, memory_id),
  shapes("shapes_generic", id),
  shapes_derivatives("shapes_derivatives_generic", id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void ShapeLagrange<_ek_igfem>::printself(std::ostream & stream, int indent) const {
  // std::string space;
  // for(Int i = 0; i < indent; i++, space += AKANTU_INDENT);

  // stream << space << "Shapes Lagrange [" << std::endl;
  // ShapeFunctions::printself(stream, indent + 1);
  // shapes.printself(stream, indent + 1);
  // shapes_derivatives.printself(stream, indent + 1);
  // stream << space << "]" << std::endl;
}

__END_AKANTU__

#endif



