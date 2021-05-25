
.PHONY: pypi help test

help:
	@echo "This is a makefile to push to pypi."
	@echo "Use make pypi to push to pypi."

test:
	py.test pyfinite --doctest-modules


dist: pyfinite setup.py
	python3 setup.py sdist

test_pypi: README.rst dist
	twine upload --repository testpypi dist/*

pypi: README.rst test dist
	twine upload --repository pypi dist/*

README.rst: README.md
	pandoc --from=markdown --to=rst --output=README.rst README.md

