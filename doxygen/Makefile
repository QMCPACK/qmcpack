all: dev-doc main

main:  dev.tag qmcpack.cfg
	doxygen qmcpack.cfg

dev-doc: dev/*.dox dev.cfg
	doxygen dev.cfg
	cp dev.tag output/

