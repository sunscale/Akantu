#!/bin/bash

./test_embedded_interface_model && diff matrix matrix_test >/dev/null 2>&1
