#Copyright 2008 Guillaume ANCIAUX (guillaume.anciaux@epfl.ch)
#
#This file is part of ParaViewHelper.
#
#ParaViewHelper is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#ParaViewHelper is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with ParaViewHelper.  If not, see <http://www.gnu.org/licenses/>.

CC=cc
CXX=g++
CXXFLAGS=-g -Wall


all: testcpp.o testc.o mesher.o
	${CXX} ${CXXFLAGS} -o testcpp testcpp.o mesher.o -L../src -lParaviewHelper -Wl,-rpath ../src -lz
	${CC} ${CXXFLAGS} -o testc testc.o mesher.o -L../src -lParaviewHelper -Wl,-rpath ../src -lz

testcpp.o: testcpp.cpp
	${CXX} ${CXXFLAGS} -c testcpp.cpp 

testc.o: testc.c
	${CC} ${CXXFLAGS} -c testc.c

mesher.o: mesher.c
	${CC} ${CXXFLAGS} -c mesher.c

clean:
	rm -f *.o