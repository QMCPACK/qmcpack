all: ug-doc dev-doc tut-doc main

main:  ug.tag dev.tag tutorials.tag qmcpack.cfg
	doxygen qmcpack.cfg

ug-doc: ug/*.dox ug.cfg
	doxygen ug.cfg
	cp ug.tag output/

dev-doc: dev/*.dox dev.cfg
	doxygen dev.cfg
	cp dev.tag output/

tut-doc: tutorials/*.dox tut.cfg
	doxygen tut.cfg
	cp tutorials.tag output/

ug-pdf: doxygen.sty 
	cp doxygen.sty tex/; cd tex; make pdf; cd -
