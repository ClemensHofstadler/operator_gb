.PHONY: install setup build clean

install: build setup
	sage -pip install --upgrade --no-index .

setup: requirements.txt
	pip install -r requirements.txt

build: 
	sage -python setup.py build_ext --inplace

clean:
	rm -rf __pycache__