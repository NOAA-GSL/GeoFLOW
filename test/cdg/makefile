#
#*********************************************************************
# Name         : GTest helper makefile:
# Date         : 2/12/2018 (DLR)
# Copyright    : 2018. Colorado State University. All rights reserved
# Description  : Helper makefile to access CMake build procedure
#*********************************************************************

all: 
	cd ../../build  ; \
	make         ; \
	make install ; \
	cd ../test/cdg
clean:  
	cd ../../build ; \
	make clean  ; \
	cd ../test/cdg
dist: 
	/bin/rm -rf ../../build/* ; \
	cd ../../build            ; \
	cmake ../src/cdg          ; \
	make                      ; \
	make install              ; \
	cd ../test/cdg
tdist: 
	/bin/rm -rf ../../build/* ; \
	cd ../../build            ; \
	cmake ../test/cdg         ; \
	make                      ; \
	make install              ; \
	cd ../test/cdg

