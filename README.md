
# Introduction

The `pyfinite` package is a python package for dealing with finite
fields and related mathematical operations. Also included is a generic
matrix package for doing matrix operations over generic fields. As an
illustration a Reed-Solomon erasure correcting code implementation is
provided using these tools.

Roughly speaking a "field" is a mathematical space where consistent
addition, subtraction, multiplication, and division operations are
defined. A "finite field" is a field where the number of elements is
finite. Perhaps the most familiar finite field is the Boolean field
where the elements are 0 and 1, addition (and subtraction) correspond
to XOR, and multiplication (and division) work as normal for 0 and 1.

More complicated finite fields are useful and interesting for
cryptography and erasure correcting codes.

# Usage

After you install via something like `pip install pyfinite`, the best way to get started is to look at the doctest examples in the following files:

  - `ffield.py`: See docstring for `FField` and `FElement` classes.
	- This shows you how to work with finite fields.
  - `genericmatrix.py`: See docstring for `GenericMatrix` class.
	- This shows you how to do matrix operations on a generic field.
  - `rs_code.py`: See docstring for `RSCode` class.
	- This shows you how to do Reed-Solomon erasure correcting codes.
  - `file_ecc.py`: See the top-level docstring for the `file_ecc` module.
    - Shows you how to encode a file into multiple pieces and decode from a subset of those pieces.
	
For example, after you install `pyfinite` and start the python
interpreter, do something like the following to see help on finite
fields:

```python
>>> from pyfinite import ffield
>>> help(ffield.FField)
```


# Apology

This code was written many years ago and hosted on an old MIT web site
before being moved to github. It is in need of some love. In
particular, it could use:

  1. Reworking to fix pep8/pylint warnings and generally better python style.
  2. More documentation.
  3. More examples.
  4. System to verify doctests in both python2 and python3.
	 - These have been manually verified but it would be nice to have a setup which can run tests on multiple versions of python.

To help or contribute please see the main project site at https://github.com/emin63/pyfinite.
