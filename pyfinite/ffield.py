
# Copyright Emin Martinian 2002.  See below for license terms.
# Version Control Info: $Id: ffield.py,v 1.5 2007/03/15 02:49:44 emin Exp $

"""
This package contains the FField class designed to perform calculations
in finite fields of characteristic two.  The following docstrings provide
detailed information on various topics:

  FField.__doc__    Describes the methods of the FField class and how
                    to use them.

  FElement.__doc__  Describes the FElement class and how to use it.

  fields_doc        Briefly describes what a finite field is and
                    establishes notation for further documentation.

  design_doc        Discusses the design of the FField class and attempts
                    to justify why certain decisions were made.

  license_doc       Describes the license and lack of warranty for
                    this code.

  testing_doc       Describes some tests to make sure the code is working
                    as well as some of the built in testing routines.

"""

import functools
from functools import reduce
import os
import os.path
import pickle
import random
import sys

# The following list of primitive polynomials are the Conway Polynomials
# from the list at
# http://www.math.rwth-aachen.de/~Frank.Luebeck/ConwayPol/cp2.html

g_primitive_polys = {}
g_primitive_polys_condensed = {
    1:   (1, 0),
    2:   (2, 1, 0),
    3:   (3, 1, 0),
    4:   (4, 1, 0),
    5:   (5, 2, 0),
    6:   (6, 4, 3, 1, 0),
    7:   (7, 1, 0),
    8:   (8, 4, 3, 2, 0),
    9:   (9, 4, 0),
    10:  (10, 6, 5, 3, 2, 1, 0),
    11:  (11, 2, 0),
    12:  (12, 7, 6, 5, 3, 1, 0),
    13:  (13, 4, 3, 1, 0),
    14:  (14, 7, 5, 3, 0),
    15:  (15, 5, 4, 2, 0),
    16:  (16, 5, 3, 2, 0),
    17:  (17, 3, 0),
    18:  (18, 12, 10, 1, 0),
    19:  (19, 5, 2, 1, 0),
    20:  (20, 10, 9, 7, 6, 5, 4, 1, 0),
    21:  (21, 6, 5, 2, 0),
    22:  (22, 12, 11, 10, 9, 8, 6, 5, 0),
    23:  (23, 5, 0),
    24:  (24, 16, 15, 14, 13, 10, 9, 7, 5, 3, 0),
    25:  (25, 8, 6, 2, 0),
    26:  (26, 14, 10, 8, 7, 6, 4, 1, 0),
    27:  (27, 12, 10, 9, 7, 5, 3, 2, 0),
    28:  (28, 13, 7, 6, 5, 2, 0),
    29:  (29, 2, 0),
    30:  (30, 17, 16, 13, 11, 7, 5, 3, 2, 1, 0),
    31:  (31, 3, 0),
    32:  (32, 15, 9, 7, 4, 3, 0),
    33:  (33, 13, 12, 11, 10, 8, 6, 3, 0),
    34:  (34, 16, 15, 12, 11, 8, 7, 6, 5, 4, 2, 1, 0),
    35:  (35, 11, 10, 7, 5, 2, 0),
    36:  (36, 23, 22, 20, 19, 17, 14, 13, 8, 6, 5, 1, 0),
    37:  (37, 5, 4, 3, 2, 1, 0),
    38:  (38, 14, 10, 9, 8, 5, 2, 1, 0),
    39:  (39, 15, 12, 11, 10, 9, 7, 6, 5, 2, 0),
    40:  (40, 23, 21, 18, 16, 15, 13, 12, 8, 5, 3, 1, 0),
    97:  (97, 6, 0),
    100: (100, 15, 0),
}

for n in g_primitive_polys_condensed.keys():
    g_primitive_polys[n] = [0]*(n+1)
    unity = 1
    for index in g_primitive_polys_condensed[n]:
        g_primitive_polys[n][index] = unity
    g_primitive_polys[n].reverse()


if sys.version_info[0] >= 3:
    def long(data):
        """Fake the `long` function since not needed after python 2"""
        return data


class FField:
    """
    The FField class implements a finite field calculator.  The
    following functions are provided:

    __init__
    add
    subtract
    multiply
    inverse
    divide
    find_degree
    multiply_without_reducing
    extended_euclid
    full_division
    show_coefficients
    show_polynomial
    get_random_element
    convert_list_to_element
    test_full_division
    test_inverse

    Most of these methods take integers or longs representing field
    elements as arguments and return integers representing the desired
    field elements as output.  See ffield.fields_doc for an explanation
    of the integer representation of field elements.

    Example of how to use the FField class:

>>> import ffield
>>> f = ffield.FField(5)  # create the field GF(2^5)
>>> a = 7  # field elements are denoted as integers from 0 to 2^5-1
>>> b = 15
>>> f.show_polynomial(a)  # show the polynomial representation of a
'x^2 + x^1 + 1'
>>> f.show_polynomial(b)
'x^3 + x^2 + x^1 + 1'
>>> c = f.multiply(a, b)  # multiply a and b modulo the field generator
>>> c
8
>>> f.show_polynomial(c)
'x^3'
>>> f.multiply(c, f.inverse(a)) == b  # verify multiplication works
1
>>> f.multiply(c, f.inverse(b)) == a  # verify multiplication works
1
>>> d = f.divide(c, b)  # since c = f.Multiply(a, b), d should give a
>>> d
7

    See documentation on the appropriate method for further details.
    """

    def __init__(self, n, gen=0, use_lut=-1):
        """
        This method constructs the field GF(2^p).  It takes one
        required argument, n = p, and two optional arguments, gen,
        representing the coefficients of the generator polynomial
        (of degree n) to use and use_lut describing whether to use
        a lookup table.  If no gen argument is provided, the
        Conway Polynomial of degree n is obtained from the table
        g_primitive_polys.

        If use_lut == 1 then a lookup table is used for
        computing finite field multiplies and divides.
        If use_lut == 0 then no lookup table is used.
        If use_lut == -1 (the default), then the code
        decides when a lookup table should be used.

        Note that you can look at the generator for the field object
        f by looking at f.generator.
        """
        self.n = n
        if gen:
            self.generator = gen
        else:
            self.generator = self.convert_list_to_element(g_primitive_polys[n])

        self.lut = None
        if use_lut == 1 or (use_lut == -1 and self.n < 10):  # use lookup table
            self.unity = 1
            self.inverse = self._inverse_for_small_field
            self.prepare_lut()
            self.multiply = self._lut_multiply
            self.divide = self._lut_divide
            self.inverse = functools.partial(self._lut_divide, 1)
        elif self.n < 15:
            self.unity = 1
            self.inverse = self._inverse_for_small_field
            self.multiply = self._multiply
            self.divide = self._divide
        else:  # Need to use longs for larger fields
            self.unity = long(1)
            self.inverse = self._inverse_for_big_field
            self.multiply = lambda a, b: self._multiply(long(a), long(b))
            self.divide = lambda a, b: self._divide(long(a), long(b))

    def prepare_lut(self):
        field_size = 1 << self.n
        lut_name = 'ffield.lut.' + repr(self.n)
        if os.path.exists(lut_name):
            fd = open(lut_name, 'rb')
            self.lut = pickle.load(fd)
            fd.close()
        else:
            self.lut = LUT()
            self.lut.mul_lut = list(range(field_size))
            self.lut.div_lut = list(range(field_size))
            self.lut.mul_lut[0] = [0] * field_size
            self.lut.div_lut[0] = ['NaN'] * field_size
            for i in range(1, field_size):
                self.lut.mul_lut[i] = [self._multiply(i, x)
                                       for x in range(field_size)]
                self.lut.div_lut[i] = [self._divide(i, x)
                                       for x in range(field_size)]
            fd = open(lut_name, 'wb')
            pickle.dump(self.lut, fd)
            fd.close()

    def _lut_multiply(self, i, j):
        return self.lut.mul_lut[i][j]

    def _lut_divide(self, i, j):
        return self.lut.div_lut[i][j]

    def add(self, x, y):
        """
        Adds two field elements and returns the result.
        """
        return x ^ y

    def subtract(self, x, y):
        """
        Subtracts the second argument from the first and returns
        the result.  In fields of characteristic two this is the same
        as the add method.
        """
        return self.add(x, y)

    def _multiply(self, f, v):
        """
        Multiplies two field elements (modulo the generator
        self.generator) and returns the result.

        See multiply_without_reducing if you don't want multiplication
        modulo self.generator.
        """
        m = self.multiply_without_reducing(f, v)
        return self.full_division(m, self.generator,
                                  self.find_degree(m), self.n)[1]

    def _inverse_for_small_field(self, f):
        """
        Computes the multiplicative inverse of its argument and
        returns the result.
        """
        return self.extended_euclid(1, f, self.generator,
                                    self.find_degree(f), self.n)[1]

    def _inverse_for_big_field(self, f):
        """
        Computes the multiplicative inverse of its argument and
        returns the result.
        """
        return self.extended_euclid(self.unity, long(f), self.generator,
                                    self.find_degree(long(f)), self.n)[1]

    def _divide(self, f, v):
        """
        divide(f,v) returns f * v^-1.
        """
        return self._multiply(f, self.inverse(v))

    def find_degree(self, v):
        """
        Find the degree of the polynomial representing the input field
        element v.  This takes O(degree(v)) operations.

        A faster version requiring only O(log(degree(v)))
        could be written using binary search...
        """
        if v:
            result = -1
            while v:
                v = v >> 1
                result = result + 1
            return result
        else:
            return 0

    def multiply_without_reducing(self, f, v):
        """
        Multiplies two field elements and does not take the result
        modulo self.generator.  You probably should not use this
        unless you know what you are doing; look at multiply instead.

        NOTE: If you are using fields larger than GF(2^15), you should
        make sure that f and v are longs not integers.
        """
        result = 0
        mask = self.unity
        i = 0
        while i <= self.n:
            if mask & v:
                result = result ^ f
            f = f << 1
            mask = mask << 1
            i = i + 1
        return result

    def extended_euclid(self, d, a, b, a_degree, b_degree):
        """
        Takes arguments (d, a, b, a_degree, b_degree) where d = gcd(a, b)
        and returns the result of the extended Euclid algorithm
        on (d, a, b).
        """
        if b == 0:
            return a, self.unity, 0
        else:
            floor_a_div_b, a_mod_b = self.full_division(a, b,
                                                        a_degree, b_degree)
            d, x, y = self.extended_euclid(d, b, a_mod_b, b_degree,
                                           self.find_degree(a_mod_b))
            return d, y, self.subtract(x, self._multiply(floor_a_div_b, y))

    def full_division(self, f, v, f_degree, v_degree):
        """
        Takes four arguments, f, v, f_degree, and v_degree where
        f_degree and v_degree are the degrees of the field elements
        f and v represented as a polynomials.
        This method returns the field elements a and b such that

            f(x) = a(x) * v(x) + b(x).

        That is, a is the divisor and b is the remainder, or in
        other words a is like floor(f/v) and b is like f modulo v.
        """
        result = 0
        i = f_degree
        mask = self.unity << i
        while i >= v_degree:
            if mask & f:
                result = result ^ (self.unity << (i - v_degree))
                f = self.subtract(f, v << (i - v_degree))
            i = i - 1
            mask = mask >> self.unity
        return result, f

    def show_coefficients(self, f):
        """
        Show coefficients of input field element represented as a
        polynomial in decreasing order.
        """
        f_degree = self.n

        result = []
        for i in range(f_degree, -1, -1):
            if (self.unity << i) & f:
                result.append(1)
            else:
                result.append(0)

        return result

    def show_polynomial(self, f):
        """
        Show input field element represented as a polynomial.
        """
        f_degree = self.find_degree(f)
        result = ''

        if f == 0:
            return '0'

        for i in range(f_degree, 0, -1):
            if (1 << i) & f:
                result = result + (' x^' + repr(i))
        if 1 & f:
            result = result + ' ' + repr(1)
        return result.strip().replace(' ', ' + ')

    def get_random_element(self, non_zero=0, max_degree=None):
        """
        Return an element from the field chosen uniformly at random
        or, if the optional argument nonZero is true, chosen uniformly
        at random from the non-zero elements, or, if the optional argument
        max_degree is provided, ensure that the result has degree less
        than max_degree.
        """
        if max_degree is None:
            max_degree = self.n
        if max_degree <= 1 and non_zero:
            return 1
        if max_degree < 31:
            return random.randint(non_zero != 0, (1 << max_degree) - 1)
        else:
            result = 0
            for i in range(0, max_degree):
                result = result ^ (random.randint(0, 1) << long(i))
            if non_zero and result == 0:
                return self.get_random_element(1)
            else:
                return result

    def convert_list_to_element(self, lst):
        """
        This method takes as input a binary list (e.g. [1, 0, 1, 1])
        and converts it to a decimal representation of a field element.
        For example, [1, 0, 1, 1] is mapped to 8 | 2 | 1 = 11.

        Note if the input list is of degree >= to the degree of the
        generator for the field, then you will have to call take the
        result modulo the generator to get a proper element in the
        field.
        """
        temp = map(lambda a, b: a << b, lst, range(len(lst) - 1, -1, -1))
        return reduce(lambda a, b: a | b, temp)

    def test_full_division(self):
        """
        Test the full_division function by generating random polynomials
        a(x) and b(x) and checking whether (c, d) == full_division(a, b)
        satisfies b*c + d == a
        """
        a = self.get_random_element(non_zero=1)
        b = self.get_random_element(non_zero=1)
        a_degree = self.find_degree(a)
        b_degree = self.find_degree(b)

        c, d = self.full_division(a, b, a_degree, b_degree)
        recon = self.add(d, self.multiply(c, b))
        assert (recon == a), ('test_full_division failed: a='
                              + repr(a) + ', b=' + repr(b) + ', c='
                              + repr(c) + ', d=' + repr(d) + ', recon=', recon)

    def test_inverse(self):
        """
        This function tests the inverse function by generating
        a random non-zero polynomials a(x) and checking if
        a * inverse(a) == 1.
        """
        a = self.get_random_element(non_zero=1)
        a_inv = self.inverse(a)
        prod = self.multiply(a, a_inv)
        assert 1 == prod, ('TestInverse failed:' + 'a=' + repr(a) + ', a_inv='
                           + repr(a_inv) + ', prod=' + repr(prod),
                           'gen=' + repr(self.generator))


class LUT:
    """
    Lookup table used to speed up some finite field operations.
    """
    pass


class FElement:
    """
    This class provides field elements which overload the
    +, -, *, %, //, / operators to be the appropriate field operation.
    Note that before creating FElement objects you must first
    create an FField object.  For example,

>>> import ffield
>>> f = ffield.FField(5)
>>> e1 = ffield.FElement(f, 7)
>>> e1
x^2 + x^1 + 1
>>> e2 = ffield.FElement(f, 19)
>>> e2
x^4 + x^1 + 1
>>> e3 = e1 + e2
>>> e3
x^4 + x^2
>>> e4 = e3 / e2
>>> e4
x^4 + x^3
>>> e4 * e2 == e3
1

    """

    def __init__(self, field, e):
        """
        The constructor takes two arguments, field, and e where
        field is an FField object and e is an integer representing
        an element in FField.

        The result is a new FElement instance.
        """
        self.f = e
        self.field = field

    def __add__(self, other):
        assert self.field == other.field
        return FElement(self.field, self.field.add(self.f, other.f))

    def __mul__(self, other):
        assert self.field == other.field
        return FElement(self.field, self.field.multiply(self.f, other.f))

    def __mod__(self, o):
        assert self.field == o.field
        return FElement(self.field,
                        self.field.full_division(
                            self.f, o.f, self.field.find_degree(self.f),
                            self.field.find_degree(o.f))[1])

    def __floordiv__(self, o):
        assert self.field == o.field
        return FElement(self.field,
                        self.field.full_division(
                            self.f, o.f, self.field.find_degree(self.f),
                            self.field.find_degree(o.f))[0])

    def __truediv__(self, other):
        assert self.field == other.field
        return FElement(self.field, self.field.divide(self.f, other.f))

    def __div__(self, *args, **kwargs):
        """syntactic sugar for calling self.__truediv__(*args, **kwargs)"""
        return self.__truediv__(*args, **kwargs)

    def __str__(self):
        return self.field.show_polynomial(self.f)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        assert self.field == other.field
        return self.f == other.f


def full_test(tests_per_field=10, size_list=None):
    """
    This function runs test_inverse and test_full_division for tests_per_field
    random field elements for each field size in size_list.  For example,
    if size_list = (1, 5, 7), then tests are run on GF(2), GF(2^5), and
    GF(2^7).  If size_list is None (which is the default), then every
    field is tested.
    """
    if size_list is None:
        size_list = g_primitive_polys.keys()
    for i in size_list:
        f = FField(i)
        for j in range(tests_per_field):
            f.test_inverse()
            f.test_full_division()


fields_doc = """
Roughly speaking a finite field is a finite collection of elements
where most of the familiar rules of math work.  Specifically, you
can add, subtract, multiply, and divide elements of a field and
continue to get elements in the field.  This is useful because
computers usually store and send information in fixed size chunks.
Thus many useful algorithms can be described as elementary operations
(e.g. addition, subtract, multiplication, and division) of these chunks.

Currently this package only deals with fields of characteristic 2.  That
is all fields we consider have exactly 2^p elements for some integer p.
We denote such fields as GF(2^p) and work with the elements represented
as p-1 degree polynomials in the indeterminate x.  That is an element of
the field GF(2^p) looks something like

     f(x) = c_{p-1} x^{p-1} + c_{p-2} x^{p-2} + ... + c_0

where the coefficients c_i are in binary.

Addition is performed by simply adding coefficients of degree i
modulo 2.  For example, if we have two field elements f and v
represented as f(x) = x^2 + 1 and v(x) = x + 1 then s = f + v
is given by (x^2 + 1) + (x + 1) = x^2 + x.  Multiplication is
performed modulo a p degree generator polynomial g(x).
For example, if f and v are as in the above example, then s = s * v
is given by (x^2 + 1) + (x + 1) mod g(x).  Subtraction turns out
to be the same as addition for fields of characteristic 2.  Division
is defined as f / v = f * v^-1 where v^-1 is the multiplicative
inverse of v.  Multiplicative inverses in groups and fields
can be calculated using the extended Euclid algorithm.

Roughly speaking the intuition for why multiplication is
performed modulo g(x), is because we want to make sure s * v
returns an element in the field.  Elements of the field are
polynomials of degree p-1, but regular multiplication could
yield terms of degree greater than p-1.  Therefore we need a
rule for 'reducing' terms of degree p or greater back down
to terms of degree at most p-1.  The 'reduction rule' is
taking things modulo g(x).

For another way to think of
taking things modulo g(x) as a 'reduction rule', imagine
g(x) = x^7 + x + 1 and we want to take some polynomial,
f(x) = x^8 + x^3 + x, modulo g(x).  We can think of g(x)
as telling us that we can replace every occurrence of
x^7 with x + 1.  Thus f(x) becomes x * x^7 + x^3 + x which
becomes x * (x + 1) + x^3 + x = x^3 + x^2 .  Essentially, taking
polynomials mod x^7 by replacing all x^7 terms with x + 1 will
force down the degree of f(x) until it is below 7 (the leading power
of g(x).  See a book on abstract algebra for more details.
"""

design_doc = """
The FField class implements a finite field calculator for fields of
characteristic two.  This uses a representation of field elements
as integers and has various methods to calculate the result of
adding, subtracting, multiplying, dividing, etc. field elements
represented AS INTEGERS OR LONGS.

The FElement class provides objects which act like a new kind of
numeric type (i.e. they overload the +, -, *, %, //, / operators,
and print themselves as polynomials instead of integers).

Use the FField class for efficient storage and calculation.
Use the FElement class if you want to play around with finite
field math the way you would in something like Matlab or
Mathematica.

--------------------------------------------------------------------
                           WHY PYTHON?

You may wonder why a finite field calculator written in Python would
be useful considering all the C/C++/Java code already written to do
the same thing (and probably faster too).  The goals of this project
are as follows, please keep them in mind if you make changes:

o  Provide an easy to understand implementation of field operations.
   Python lends itself well to comments and documentation.  Hence,
   we hope that in addition to being useful by itself, this project
   will make it easier for people to implement finite field
   computations in other languages.  If you've ever looked at some
   of the highly optimized finite field code written in C, you will
   understand the need for a clear reference implementation of such
   operations.

o  Provide easy access to a finite field calculator.
   Since you can just start up the Python interpreter and do
   computations, a finite field calculator in Python lets you try
   things out, check your work for other algorithms, etc.
   Furthermore since a wealth of numerical packages exist for python,
   you can easily write simulations or algorithms which draw upon
   such routines with finite fields.

o  Provide a platform independent framework for coding in Python.
   Many useful error control codes can be implemented based on
   finite fields.  Some examples include error/erasure correction,
   cyclic redundancy checks (CRCs), and secret sharing.  Since
   Python has a number of other useful Internet features being able
   to implement these kinds of codes makes Python a better framework
   for network programming.

o  Leverages Python arbitrary precision code for large fields.
   If you want to do computations over very large fields, for example
   GF(2^p) with p > 31 you have to write lots of ugly bit field
   code in most languages.  Since Python has built in support for
   arbitrary precision integers, you can make this code work for
   arbitrary field sizes provided you operate on longs instead of
   ints.  In python 3 this is automatic, but in python2 you need
   to give as input numbers like 0L, 1L, 1L << 55, etc. to force
   using longs.

--------------------------------------------------------------------
                            BASIC DESIGN


The basic idea is to index entries in the finite field of interest
using integers and design the class methods to work properly on this
representation.  Using integers is efficient since integers are easy
to store and manipulate and allows us to handle arbitrary field sizes
without changing the code if we instead switch to using longs.

Specifically, an integer represents a bit string

  c = c_{p-1} c_{p-2} ... c_0.

which we interpret as the coefficients of a polynomial representing a
field element

  f(x) = c_{p-1} x^{p-1} + c_{p-2} x^{p-2} + ... + c_0.

--------------------------------------------------------------------
                             FUTURE
In the future, support for fields of other
characteristic may be added (if people want them).  Since computers
have built in parallelized operations for fields of characteristic
two (i.e. bitwise and, or, xor, etc.), this implementation uses
such operations to make most of the computations efficient.

"""


license_doc = """
  This code was originally written by Emin Martinian (emin@allegro.mit.edu).
  You may copy, modify, redistribute in source or binary form as long
  as credit is given to the original author.  Specifically, please
  include some kind of comment or docstring saying that Emin Martinian
  was one of the original authors.  Also, if you publish anything based
  on this work, it would be nice to cite the original author and any
  other contributors.

  There is NO WARRANTY for this software just as there is no warranty
  for GNU software (although this is not GNU software).  Specifically
  we adopt the same policy towards warranties as the GNU project:

  BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE PROGRAM 'AS IS' WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED
OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS
TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE
PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
REPAIR OR CORRECTION.

  IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES,
INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING
OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED
TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY
YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER
PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE
POSSIBILITY OF SUCH DAMAGES.
"""

testing_doc = """
The FField class has a number of built in testing functions such as
test_full_division, test_inverse.  The simplest thing to
do is to call the full_test method.

>>> import ffield
>>> ffield.full_test(size_list=None, tests_per_field=100)

# To decrease the testing time you can either decrease the tests_per_field
# or you can only test the field sizes you care about by doing something
# like size_list = [2, 7, 20] in the ffield.full_test command above.

If any problems occur, assertion errors are raised.  Otherwise
nothing is returned.  Note that you can also use the doctest
package to test all the python examples in the documentation
by typing 'python ffield.py' or 'python -v ffield.py' at the
command line.
"""


# The following code is used to make the doctest package
# check examples in docstrings.
__test__ = {
    'testing_doc': testing_doc
}


def _test():
    import doctest
    import ffield

    return doctest.testmod(ffield)


if __name__ == "__main__":
    print('Starting automated tests (this may take a while)')
    _test()
    print('Tests finished.')
