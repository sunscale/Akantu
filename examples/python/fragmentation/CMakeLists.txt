add_mesh(fragmentation_mesh fragmentation_mesh.geo DIM 2 ORDER 1)

register_example(fragmentation
  SCRIPT fragmentation.py
  FILES_TO_COPY material.dat
  DEPENDS fragmentation_mesh
  PYTHON
  )
