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
clear:  
	cd ../../build ; \
	/bin/rm -r  *; \
	cd ../test/cdg
docs:
	/bin/rm -rf ../../build/* ; \
	cd ../../build            ; \
	cmake -D GCONF=duane        \
              ../                 ; \
	make docs
dist: 
	/bin/rm -rf ../../build/* ; \
	cd ../../build            ; \
	cmake -D GCONF=duane        \
              ../                 ; \
	make                      ; \
	make install              ; \
	cd ../test/cdg

