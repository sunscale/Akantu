#!/bin/bash

svn info | grep "vision" | head -n 1 | cut -d ":" -f 2