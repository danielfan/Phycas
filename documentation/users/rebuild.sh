#!/bin/sh
set -x 
python regenerate_inputs.py 
pdflatex manual.tex 
bibtex manual
pdflatex manual.tex 
pdflatex manual.tex 
pdflatex manual.tex 
