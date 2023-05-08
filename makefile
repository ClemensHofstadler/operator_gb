.PHONY: install build clean

install: build
	sage -pip install --upgrade .

build: 
	sage -python setup.py build_ext --inplace

clean:
	rm -rf __pycache__