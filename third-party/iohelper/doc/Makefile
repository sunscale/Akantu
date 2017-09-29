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


manual:
	latex manual.tex
	latex manual.tex
	dvips -Ppdf -G0 manual.dvi -o manual.ps
#	latex2html -local_icons manual.tex 
	ps2pdf -dPDFsettings=/prepress manual.ps manual.pdf


clean:
	rm -f *.aux *.log *.pdf *.ps *.dvi *.bbl *.blg *.toc *.out* *~	
