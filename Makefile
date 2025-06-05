.PHONY: clean build

build:
	chmod u+x build.sh
	cotocoa_mode=COTOCOA_DISABLED ./build.sh

requester:
	chmod u+x build.sh
	cotocoa_mode=COTOCOA_REQUESTER ./build.sh

worker:
	chmod u+x build.sh
	cotocoa_mode=COTOCOA_WORKER ./build.sh

clean:
	-make -C lib/mtarm clean
	-make -C lib/ohhelp/ohhelp-1.1.1 clean
	-make -C lib/psort clean
	-make -C lib/cotocoa/src -f Makefile.env clean
	-rm -rf build/
	-rm -rf bin/
	-rm -rf *.mod
