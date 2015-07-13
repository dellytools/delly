PIP=../venv/bin/pip

all:
	$(PIP) freeze > dev-python2.7.txt
	grep -v -e ^ipython -e ^pep8 dev-python2.7.txt > python2.7.txt 

