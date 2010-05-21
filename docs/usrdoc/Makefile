default:	retypeset

bib:	qmcpack-u.bib
	bibtex qmcpack-u

clean:
	rm -f *.aux *.toc *.blg *.bbl *.log qmcpack-u.pdf

retypeset:	qmcpack-u.tex qmcpack-u-ov.tex qmcpack-u-gs.tex qmcpack-u-cnf.tex qmcpack-u-run.tex qmcpack-u-ex.tex qmcpack-u.bib
	pdflatex qmcpack-u.tex

typeset:
	make clean; pdflatex qmcpack-u.tex; make bib; pdflatex qmcpack-u.tex; pdflatex qmcpack-u.tex
