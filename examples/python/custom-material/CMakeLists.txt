register_example(bi-material
  SCRIPT bi-material.py
  FILES_TO_COPY material.dat square.geo)

register_example(custom-material
  SCRIPT custom-material.py
  FILES_TO_COPY material.dat bar.geo)
