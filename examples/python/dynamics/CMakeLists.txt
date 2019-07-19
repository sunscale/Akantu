add_mesh(bar_dynamics_mesh bar.geo DIM 2)
register_example(dynamics-python
  SCRIPT dynamics.py
  PYTHON
  FILES_TO_COPY material.dat
  DEPENDS bar_dynamics_mesh)
