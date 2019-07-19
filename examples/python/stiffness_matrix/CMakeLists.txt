add_mesh(plate plate.geo DIM 2)

register_example(stiffness_matrix
  SCRIPT stiffness_matrix.py
  FILES_TO_COPY material.dat
  DEPENDS plate
  )
