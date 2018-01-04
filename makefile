
.PHONY: pypi help

help:
	@echo "This is a makefile to push to pypi."
	@echo "Use make pypi to push to pypi."

pypi: README.rst
	 python3 setup.py sdist upload -r pypi

README.rst: README.md
	pandoc --from=markdown --to=rst --output=README.rst README.md

