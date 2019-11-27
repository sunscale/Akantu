add_mesh(plate_hole_mesh plate.geo DIM 2 ORDER 1)

register_example(plate-hole
  SCRIPT plate.py
  FILES_TO_COPY material.dat
  DEPENDS plate_hole_mesh
  PYTHON
  )
