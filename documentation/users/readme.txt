Run regenerate_inputs.py before compiling manual.tex. The regenerate_inputs.py
script creates several latex files input by manual.tex using information in
the Phycas.__init__ function in Phycas.py.

The following script can be used to automate the build process:

rm -f *.pdf *.log *.aux *.dvi *.idx *.bbl *.blg *.ilg *.ind *.out *.toc
pdflatex manual
bibtex manual
makeindex manual.idx
pdflatex manual
pdflatex manual
