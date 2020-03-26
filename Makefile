## install the python code to the system python packages

install:
	sage -python setup.py install

## for development work use sage -python setup.py develop
install-develop:
	sage -python setup.py develop
