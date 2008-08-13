#!/bin/sh
set -x 
pdflatex phycas_woods_hole_08.tex 
bibtex phycas_woods_hole_08
pdflatex phycas_woods_hole_08.tex 
pdflatex phycas_woods_hole_08.tex 
pdflatex phycas_woods_hole_08.tex 
