#!/bin/sh
set -x 
pdflatex phycas_woods_hole_09.tex 
bibtex phycas_woods_hole_09
pdflatex phycas_woods_hole_09.tex 
pdflatex phycas_woods_hole_09.tex 
pdflatex phycas_woods_hole_09.tex 
