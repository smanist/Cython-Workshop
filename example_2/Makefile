cy:
	cd ./beam/ && python setup_cy.py build_ext --inplace
	mv ./beam/beam/cyUtils.so ./beam
	rm -r ./beam/beam

pc:
	cd ./beam/ && python setup_pc.py build_ext --inplace

so:
	cd ./beam/ && python setup_so.py build_ext --inplace
	cd ./beam/ && gcc -shared -O3 -Wall -fPIC -c pc_utils.c -o pc_utils.so
	#gcc -shared -O3 -Wall -fPIC -fopenmp -c pc_utils.c -o pc_utils.so

clean:
	rm ./beam/pcUtils.c || continue
	rm ./beam/cyUtils.c || continue
	rm ./beam/*.so || continue
	rm ./beam/*.pyc || continue
	rm -r ./beam/build || continue

