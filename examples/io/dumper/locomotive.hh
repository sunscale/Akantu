using namespace akantu;

/* -------------------------------------------------------------------------- */
void applyRotation(Array<Real> center, Real angle, Array<Real> & nodes,
		   Array<Real> * displacement, Array<UInt> node_group) {

  for (UInt i = 0; i < node_group.getSize(); ++i) {
    
    Array<Real> pos_rel(center.getSize());

    for (UInt j = 0; j < center.getSize(); ++j) {
      
      pos_rel(j) = nodes(node_group(i),j) - center(j);
    }

    Real radius = std::sqrt(pos_rel[0]*pos_rel[0]+pos_rel[1]*pos_rel[1]);

    if (std::abs(radius) < Math::getTolerance()) continue;

    Real phi_i = std::acos(pos_rel[0]/radius);

    if (pos_rel[1] < 0) phi_i *= -1;

    (*displacement)(node_group(i),0) = std::cos(phi_i-angle)*radius - pos_rel[0];
    (*displacement)(node_group(i),1) = std::sin(phi_i-angle)*radius - pos_rel[1];
  }
}
