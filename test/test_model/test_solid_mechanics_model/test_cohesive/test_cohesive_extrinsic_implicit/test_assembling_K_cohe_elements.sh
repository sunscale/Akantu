#!/bin/bash

./test_assembling_K_cohe_elements && diff K_matrix_verified K_matrix_test >/dev/null 2>&1
