PIP=../venv/bin/pip

all:
	$(PIP) freeze | grep numpy > python2.7.txt
	$(PIP) freeze | grep -v numpy >> python2.7.txt
