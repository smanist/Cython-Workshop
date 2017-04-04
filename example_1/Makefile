all:
	cd ./prime/ && python setup_eg1.py build_ext --inplace
	mv ./prime/prime/prime_cython.so ./prime/
	rm -r ./prime/prime/

clean:
	rm ./prime/*.c || continue
	rm ./prime/*.so || continue
	rm ./prime/*.pyc || continue
	rm -r ./prime/build || continue
