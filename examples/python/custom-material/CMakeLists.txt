add_mesh(square square.geo DIM 2)
register_example(bi-material
  SCRIPT bi-material.py
  PYTHON
  FILES_TO_COPY material.dat
  DEPENDS square)

add_mesh(bar bar.geo DIM 2)
register_example(custom-material
  SCRIPT custom-material.py
  PYTHON
  FILES_TO_COPY material.dat
  DEPENDS bar)
