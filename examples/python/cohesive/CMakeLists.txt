add_mesh(plate plate.geo DIM 2)

register_example(cohesive_python
  SCRIPT cohesive.py
  PYTHON
  FILES_TO_COPY material.dat
  DEPENDS plate
  )
