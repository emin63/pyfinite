Introduction
============

The ``pyfinite`` package is a python package for dealing with finite
fields and related mathematical operations. Also included is a generic
matrix package for doing matrix operations over generic fields. As an
illustration a Reed-Solomon erasure correcting code implementation is
provided using these tools.

Roughly speaking a "field" is a mathematical space where consistent
addition, subtraction, multiplication, and division operations are
defined. A "finite field" is a field where the number of elements is
finite. Perhaps the most familiar finite field is the Boolean field
where the elements are 0 and 1, addition (and subtraction) correspond to
XOR, and multiplication (and division) work as normal for 0 and 1.

More complicated finite fields are useful and interesting for
cryptography and erasure correcting codes.

Usage
=====

After you install via something like ``pip install pyfinite``, the best
way to get started is to look at the doctest examples in the following
files:

-  ``ffield.py``: See docstring for ``FField`` and ``FElement`` classes.

   -  This shows you how to work with finite fields.

-  ``genericmatrix.py``: See docstring for ``GenericMatrix`` class.

   -  This shows you how to do matrix operations on a generic field.

-  ``rs_code.py``: See docstring for ``RSCode`` class.

   -  This shows you how to do Reed-Solomon erasure correcting codes.

-  ``file_ecc.py``: See the top-level docstring for the ``file_ecc``
   module.

   -  Shows you how to encode a file into multiple pieces and decode
      from a subset of those pieces.

For example, after you install ``pyfinite`` and start the python
interpreter, do something like the following to see help on finite
fields:

.. code:: python

    >>> from pyfinite import ffield
    >>> help(ffield.FField)

or if you want to dive right in, you can try something like the
following:

.. code:: python

    >>> from pyfinite import ffield
    >>> F = ffield.FField(5) # create the field GF(2^5)
    >>> a = 7    # field elements are denoted as integers from 0 to 2^5-1
    >>> b = 15
    >>> F.ShowPolynomial(a) # show the polynomial representation of a
    'x^2 + x^1 + 1'
    >>> c = F.Multiply(a,b) # multiply a and b modulo the field generator
    >>> c
    8
    >>> F.ShowPolynomial(c)
    'x^3'

Alternatively, you can jump into the ``genericmatrix.py`` package with
something like:

.. code:: python

    >>> import genericmatrix
    >>> v = genericmatrix.GenericMatrix((3,3))
    >>> v.SetRow(0,[0.0, -1.0, 1.0])
    >>> v.SetRow(1,[1.0, 1.0, 1.0])
    >>> v.SetRow(2,[1.0, 1.0, -1.0])
    >>> v
    <matrix
      0.0 -1.0  1.0
      1.0  1.0  1.0
      1.0  1.0 -1.0>
    >>> vi = v.Inverse()

Then for some real fun, you can try experimenting with generic matrix
operations on elements of a finite field! The nice thing about the
``genericmatrix`` module is that it only relies on the standard python
arithmetic operators so you can use it for anything with sane ``+``,
``-``, ``*``, and ``/`` operators. See the help on ``genericmatrix`` for
more info.

Finally, if you just want erasure correction, see the docs for the
``rs_code`` and ``file_ecc`` modules via something like

.. code:: python

    >>> import rs_code, file_ecc
    >>> help(file_ecc)
    >>> help(rs_code)

Future work
===========

This code was written many years ago and hosted on an old MIT web site
under the name ``py_ecc`` before being moved to github. It is in need of
some love. In particular, it could use:

1. Reworking to fix pep8/pylint warnings and generally better python
   style.
2. More documentation.
3. More examples.
4. Travis setup to verify doctests in both python2 and python3.

   -  These have been manually verified but it would be nice to have a
      setup which can run tests on multiple versions of python in an
      automated way.

To help or contribute please see the main project site at
https://github.com/emin63/pyfinite.
