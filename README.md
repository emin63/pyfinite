
# Introduction

The `pyfinite` package is a python package for dealing with finite
fields and related mathematical operations. Also included is a generic
matrix package for doing matrix operations over generic fields. As an
illustraction a Reed-Solomon erasure correcting code implementation is
provided using these tools.

Roughly speaking a "field" is a mathematical space where consistent
addition, subtraction, multiplication, and division operations are
defined. A "finite field" is a field where the number of elements is
finite. Perhaps the most familiar finite field is the boolean field
where the elements are 0 and 1, addition (and subtraction) correspond
to XOR, and multiplication (and division) work as normal for 0 and 1.

More complicated finite fields are useful and interesting for
cryptography and erasure correcting codes.

# Apology

This code was written many years ago and hosted on an old MIT web site
before being moved to github. It is in need of some love. In
particular, it could use:

  1. Reworking to fix pep8/pylint warninsg and generally better python style.
  2. A setup.py file so that people can install via pip or setuptools.
  3. More documentation.
  4. More examples.

To help or contribute please see the main project site at https://github.com/emin63/pyfinite.
