
# Copyright Emin Martinian 2002.  See below for license terms.
# Version Control Info: $Id: genericmatrix.py,v 1.5 2008-01-05 22:08:44 emin Exp $

"""
This package implements the GenericMatrix class to provide matrix
operations for any type that supports the multiply, add, subtract,
and divide operators.  For example, this package can be used to
do matrix calculations over finite fields using the ffield package
available at http://martinian.com.

The following docstrings provide detailed information on various topics:

  GenericMatrix.__doc__   Describes the methods of the GenericMatrix
                          class and how to use them.

  license_doc             Describes the license and lack of warranty
                          for this code.

  testing_doc             Describes some tests to make sure the code works.

"""

from functools import reduce
from itertools import chain
import operator


class GenericMatrix:

    """
    The GenericMatrix class implements a matrix with works with
    any generic type supporting addition, subtraction, multiplication,
    and division.  Matrix multiplication, addition, and subtraction
    are implemented as are methods for finding inverses,
    LU (actually LUP) decompositions, and determinants.  A complete
    list of user callable methods is:

    __init__
    __repr__
    __mul__
    __add__
    __sub__
    __setitem__
    __getitem__
    size
    set_row
    get_row
    get_column
    copy
    make_similar_matrix
    swap_rows
    mul_row
    add_row
    add_col
    mul_add_row
    left_mul_column_vec
    lower_gaussian_elim
    inverse
    determinant
    lup

    A quick and dirty example of how to use the GenericMatrix class
    for matricies of floats is provided below.
    
>>> import genericmatrix
>>> v = genericmatrix.GenericMatrix((3, 3))
>>> v.set_row(0, [0.0, -1.0,  1.0])
>>> v.set_row(1, [1.0,  1.0,  1.0])
>>> v.set_row(2, [1.0,  1.0, -1.0])
>>> v
<matrix
  0.0 -1.0  1.0
  1.0  1.0  1.0
  1.0  1.0 -1.0>
>>> vi = v.inverse()
>>> vi
<matrix
  1.0  0.0  1.0
 -1.0  0.5 -0.5
 -0.0  0.5 -0.5>
>>> (vi * v) - v.make_similar_matrix(v.size(), 'i')
<matrix
 0.0 0.0 0.0
 0.0 0.0 0.0
 0.0 0.0 0.0>

# See what happens when we try to invert a non-invertible matrix

>>> v[0, 1] = 0.0
>>> v
<matrix
  0.0  0.0  1.0
  1.0  1.0  1.0
  1.0  1.0 -1.0>
>>> abs(v.determinant())
0.0
>>> v.inverse()
Traceback (most recent call last):
     ...
ValueError: matrix not invertible

# LUP decomposition will still work even if inverse() won't.

>>> l, u, p = v.lup()
>>> l
<matrix
 1.0 0.0 0.0
 0.0 1.0 0.0
 1.0 0.0 1.0>
>>> u
<matrix
  1.0  1.0  1.0
  0.0  0.0  1.0
  0.0  0.0 -2.0>
>>> p
<matrix
 0.0 1.0 0.0
 1.0 0.0 0.0
 0.0 0.0 1.0>
>>> p * v - l * u
<matrix
 0.0 0.0 0.0
 0.0 0.0 0.0
 0.0 0.0 0.0>

# Operate on some column vectors using v.
# The LeftMulColumnVec methods lets us do this without having
# to construct a new GenericMatrix to represent each column vector.
>>> v.left_mul_column_vec([1.0,  2.0, 3.0])
[3.0, 6.0, 0.0]
>>> v.left_mul_column_vec([1.0, -2.0, 1.0])
[1.0, 0.0, -2.0]

# Most of the stuff above could be done with something like matlab.
# But, with this package you can do matrix ops for finite fields.
>>> XOR = lambda x, y: x^y
>>> AND = lambda x, y: x&y
>>> DIV = lambda x, y: x
>>> m = GenericMatrix(size=(3, 4), zero_element=0, identity_element=1, add=XOR, mul=AND, sub=XOR, div=DIV)
>>> m.set_row(0, [0, 1, 0, 0])
>>> m.set_row(1, [0, 1, 0, 1])
>>> m.set_row(2, [0, 0, 1, 0])
>>> # You can't invert m since it isn't square, but you can still 
>>> # get the LUP decomposition or solve a system of equations.
>>> l, u, p = v.lup()
>>> p*v-l*u
<matrix
 0.0 0.0 0.0
 0.0 0.0 0.0
 0.0 0.0 0.0>
>>> b = [1, 0, 1]
>>> x = m.solve(b)
>>> b == m.left_mul_column_vec(x)
1

    """

    def __init__(self, size=(2, 2), zero_element=0.0, identity_element=1.0,
                 add=operator.__add__, sub=operator.__sub__,
                 mul=operator.__mul__, div=operator.__truediv__,
                 eq=operator.__eq__, f_str=lambda x: repr(x),
                 equals_zero=None, fill_mode='z'):
        """
        Function:     __init__(size, zero_element, identity_element, add, sub,
                               mul, div, eq, f_str, equals_zero, fill_mode)

        Description:  This is the constructor for the GenericMatrix
                      class.  All arguments are optional and default
                      to producing a 2-by-2 zero matrix for floats.
                      A detailed description of arguments follows:

             size: A tuple of the form (numRows, numColumns)

             zero_element: An object representing the additive
                           identity (i.e. 'zero') for the data
                           type of interest.
                   
             identity_element: An object representing the multiplicative
                               identity (i.e. 'one') for the data
                               type of interest.

             add, sub, mul, div: Functions implementing basic arithmetic
                                 operations for the type of interest.

             eq: A function such that eq(x, y) == 1 if and only if x == y.

             f_str: A function used to produce a string representation of
                    the type of interest.

             equals_zero: A function used to decide if an element is
                          essentially zero.  For floats, you could use
                          lambda x: abs(x) < 1e-6.

             fill_mode: This can either be 'e' in which case the contents
                        of the matrix are left empty, 'z', in which case
                        the matrix is filled with zeros, 'i' in which
                        case an identity matrix is created, or a two
                        argument function which is called with the row
                        and column of each index and produces the value
                        for that entry.  Default is 'z'.
        """
        if equals_zero is None:
            def equals_zero(x):
                return self.eq(self.zero_element, x)

        self.equals_zero = equals_zero
        self.add = add
        self.sub = sub
        self.mul = mul
        self.div = div
        self.eq = eq
        self.str = f_str
        self.zero_element = zero_element
        self.identity_element = identity_element
        self.rows, self.cols = size
        self.data = []

        if fill_mode == 'e':
            return
        elif fill_mode == 'z':
            def fill_mode(_1, _2):
                return self.zero_element
        elif fill_mode == 'i':
            def fill_mode(x, y):
                return self.identity_element if self.eq(x, y) \
                    else self.zero_element

        for i in range(self.rows):
            self.data.append(list(map(
                fill_mode, [i] * self.cols, range(self.cols))))

    def make_similar_matrix(self, size, fill_mode):
        """
        make_similar_matrix(self, size, fill_mode)

        Return a matrix of the given size filled according to fill_mode
        with the same zero_element, identity_element, add, sub, etc.
        as self.

        For example, self.make_similar_matrix(self.size(), 'i') returns
        an identity matrix of the same shape as self.
        """
        return GenericMatrix(size=size, zero_element=self.zero_element,
                             identity_element=self.identity_element,
                             add=self.add, sub=self.sub,
                             mul=self.mul, div=self.div, eq=self.eq,
                             f_str=self.str, equals_zero=self.equals_zero,
                             fill_mode=fill_mode)

    def __repr__(self):
        m = 0
        # find the fattest element
        for r in self.data:
            for c in r:
                le = len(self.str(c))
                if le > m:
                    m = le
        f = '%%%ds' % (m+1)
        s = '<matrix'
        for r in self.data:
            s = s + '\n'
            for c in r:
                s = s + (f % self.str(c))
        s = s + '>'
        return s

    def __mul__(self, other):
        if self.cols != other.rows:
            raise ValueError("dimension mismatch")

        result = self.make_similar_matrix((self.rows, other.cols), 'z')
        for i in range(self.rows):
            for j in range(other.cols):
                result.data[i][j] = reduce(self.add, map(
                    self.mul, self.data[i], other.get_column(j)))
        return result

    def __add__(self, other):
        if self.cols != other.cols or self.rows != other.rows:
            raise ValueError("dimension mismatch")
        result = self.make_similar_matrix(size=self.size(), fill_mode='z')
        for i in range(self.rows):
            for j in range(other.cols):
                result.data[i][j] = self.add(self.data[i][j], other.data[i][j])
        return result
    
    def __sub__(self, other):
        if self.cols != other.cols or self.rows != other.rows:
            raise ValueError("dimension mismatch")
        result = self.make_similar_matrix(size=self.size(), fill_mode='z')
        for i in range(self.rows):
            for j in range(other.cols):
                result.data[i][j] = self.sub(self.data[i][j],
                                             other.data[i][j])
        return result

    def __setitem__(self, coords, data):
        """__setitem__((x, y), data) sets item row x and column y to data."""
        x_coord, y_coord = coords
        self.data[x_coord][y_coord] = data

    def __getitem__(self, coords):
        """__getitem__((x, y)) gets item at row x and column y."""
        x_coord, y_coord = coords
        return self.data[x_coord][y_coord]

    def size(self):
        """returns (rows, columns)"""
        return len(self.data), len(self.data[0])

    def set_row(self, r, result):
        """set_row(r, result) sets row r to result."""
        assert len(result) == self.cols, ('Wrong # columns in row: ' +
                                          'expected ' + repr(self.cols) +
                                          ', got ' + repr(len(result)))
        self.data[r] = list(result)

    def get_row(self, r):
        """get_row(r) returns a copy of row r."""
        return list(self.data[r])

    def get_column(self, c):
        """get_column(c) returns a copy of column c."""
        if c >= self.cols:
            raise ValueError('matrix does not have that many columns')
        result = []
        for r in self.data:
            result.append(r[c])
        return result

    def transpose(self):
        old_data = self.data
        self.data = []
        for r in range(self.cols):
            self.data.append([])
            for c in range(self.rows):
                self.data[r].append(old_data[c][r])
        rows = self.rows
        self.rows = self.cols
        self.cols = rows

    def copy(self):
        result = self.make_similar_matrix(size=self.size(), fill_mode='e')

        for r in self.data:
            result.data.append(list(r))
        return result

    def sub_matrix(self, row_start, row_end, col_start=0, col_end=None):
        """
        sub_matrix(self, row_start, row_end, col_start, col_end)
        Create and return a sub matrix containing rows
        row_start through row_end (inclusive) and columns
        col_start through col_end (inclusive).
        """
        if not col_end:
            col_end = self.cols - 1
        if row_end >= self.rows:
            raise ValueError('row_end too big: row_end >= self.rows')
        result = self.make_similar_matrix((row_end - row_start + 1,
                                           col_end - col_start + 1), 'e')

        for i in range(row_start, row_end + 1):
            result.data.append(list(self.data[i][col_start:(col_end + 1)]))

        return result
 
    def un_sub_matrix(self, row_start, row_end, col_start, col_end):
        """
        un_sub_matrix(self, row_start, row_end, col_start, col_end)
        Create and return a sub matrix containing everything except
        rows row_start through row_end (inclusive)
        and columns col_start through col_end (inclusive).
        """
        result = self.make_similar_matrix((self.rows-(row_end-row_start),
                                           self.cols-(col_end-col_start)), 'e')

        for i in chain(range(0, row_start), range(row_end, self.rows)):
            result.data.append(list(self.data[i][0:col_start] +
                                    self.data[i][col_end:]))

        return result

    def swap_rows(self, i, j):
        temp = list(self.data[i])
        self.data[i] = list(self.data[j])
        self.data[j] = temp

    def mul_row(self, r, m, start=0):
        """
        Function: mul_row(r, m, start=0)
        Multiply row r by m starting at optional column start (default 0).
        """
        row = self.data[r]
        for i in range(start, self.cols):
            row[i] = self.mul(row[i], m)

    def add_row(self, i, j):
        """
        Add row i to row j.
        """
        self.data[j] = list(map(self.add, self.data[i], self.data[j]))

    def add_col(self, i, j):
        """
        Add column i to column j.
        """
        for r in range(self.rows):
            self.data[r][j] = self.add(self.data[r][i], self.data[r][j])

    def mul_add_row(self, m, i, j):
        """
        Multiply row i by m and add to row j.
        """
        self.data[j] = list(map(self.add,
                                map(self.mul, [m]*self.cols, self.data[i]),
                                self.data[j]))

    def left_mul_column_vec(self, col_vec):
        """
        Function:       left_mul_column_vec(c)
        Purpose:        Compute the result of self * c.
        Description:    This function takes as input a list c,
                        computes the desired result and returns it
                        as a list.  This is sometimes more convenient
                        than constructed a new GenericMatrix to represent
                        c, computing the result and extracting c to a list.
        """
        if self.cols != len(col_vec):
            raise ValueError('dimension mismatch')
        result = list(range(self.rows))
        for r in range(self.rows):
            result[r] = reduce(self.add, map(self.mul, self.data[r], col_vec))
        return result

    def find_row_leader(self, start_row, c):
        for r in range(start_row, self.rows):
            if not self.eq(self.zero_element, self.data[r][c]):
                return r
        return -1

    def find_col_leader(self, r, start_col):
        for c in range(start_col, self.cols):
            if not self.equals_zero(self.data[r][c]):
                return c
        return -1    

    def partial_lower_gauss_elim(self, row_index, col_index, result_inv):
        """
        Function: partial_lower_gauss_elim(row_index, col_index, result_inv)
        
        This function does partial Gaussian elimination on the part of
        the matrix on and below the main diagonal starting from
        row_index.  In addition to modifying self, this function
        applies the required elementary row operations to the input
        matrix result_inv.

        By partial, what we mean is that if this function encounters
        an element on the diagonal which is 0, it stops and returns
        the corresponding row_index.  The caller can then permute
        self or apply some other operation to eliminate the zero
        and recall partial_lower_gauss_elim.
        
        This function is meant to be combined with upper_inverse
        to compute inverses and LU decompositions.
        """
        last_row = self.rows-1
        while row_index < last_row:
            if col_index >= self.cols:
                return row_index, col_index
            if self.eq(self.zero_element, self.data[row_index][col_index]):
                # self[row_index,col_index] = 0 so quit.
                return row_index, col_index
            divisor = self.div(self.identity_element,
                               self.data[row_index][col_index])
            for k in range(row_index + 1, self.rows):
                next_term = self.data[k][col_index]
                if self.zero_element != next_term:
                    multiple = self.mul(divisor, self.sub(self.zero_element,
                                                          next_term))
                    self.mul_add_row(multiple, row_index, k)
                    result_inv.mul_add_row(multiple, row_index, k)
            row_index = row_index + 1
            col_index = col_index + 1
        return row_index, col_index

    def lower_gaussian_elim(self, result_inv=''):
        """
        Function:       lower_gaussian_elim(r)
        Purpose:        Perform Gaussian elimination on self to eliminate
                        all terms below the diagonal.
        Description:    This method modifies self via Gaussian elimination
                        and applies the elementary row operations used in
                        this transformation to the input matrix, r
                        (if one is provided, otherwise a matrix with
                        identity elements on the main diagonal is
                        created to serve the role of r).

                        Thus if the input, r, is an identity matrix, after
                        the call it will represent the transformation
                        made to perform Gaussian elimination.

                        The matrix r is returned.
        """
        if result_inv == '':
            result_inv = self.make_similar_matrix(self.size(), 'i')

        row_index, col_index = 0, 0
        last_row = min(self.rows - 1, self.cols)
        last_col = self.cols - 1
        while row_index < last_row and col_index < last_col:
            leader = self.find_row_leader(row_index, col_index)
            if leader < 0:
                col_index = col_index + 1
                continue
            if leader != row_index:
                result_inv.add_row(leader, row_index)
                self.add_row(leader, row_index)
            row_index, col_index = self.partial_lower_gauss_elim(
                row_index, col_index, result_inv)
        return result_inv

    def upper_inverse(self, result_inv=''):
        """
        Function: upper_inverse(result_inv)
        
        Assumes that self is an upper triangular matrix like

          [a b c ... ]
          [0 d e ... ]
          [0 0 f ... ]
          [.     .   ]
          [.      .  ]
          [.       . ]

        and performs Gaussian elimination to transform self into
        the identity matrix.  The required elementary row operations
        are applied to the matrix result_inv passed as input.  For
        example, if the identity matrix is passed as input, then the
        value returned is the inverse of self before the function
        was called.

        If no matrix, result_inv, is provided as input then one is
        created with identity elements along the main diagonal.
        In either case, result_inv is returned as output.
        """
        if result_inv == '':
            result_inv = self.make_similar_matrix(self.size(), 'i')
        last_col = min(self.rows, self.cols)
        for colIndex in range(0, last_col):
            if self.zero_element == self.data[colIndex][colIndex]:
                raise ValueError('matrix not invertible')
            divisor = self.div(self.identity_element,
                               self.data[colIndex][colIndex])
            if self.identity_element != divisor:
                self.mul_row(colIndex, divisor, colIndex)
                result_inv.mul_row(colIndex, divisor)
            for row_to_elim in range(0, colIndex):
                multiple = self.sub(self.zero_element,
                                    self.data[row_to_elim][colIndex])
                self.mul_add_row(multiple, colIndex, row_to_elim)
                result_inv.mul_add_row(multiple, colIndex, row_to_elim)
        return result_inv

    def inverse(self):
        """
        Function:       inverse
        Description:    Returns the inverse of self without modifying
                        self.  An exception is raised if the matrix
                        is not invertible.
        """

        working_copy = self.copy()
        result = self.make_similar_matrix(self.size(), 'i')
        working_copy.lower_gaussian_elim(result)
        working_copy.upper_inverse(result)
        return result

    def determinant(self):
        """
        Function:       determinant
        Description:    Returns the determinant of the matrix or raises
                        a ValueError if the matrix is not square.
        """
        if self.rows != self.cols:
            raise ValueError('matrix not square')
        working_copy = self.copy()
        result = self.make_similar_matrix(self.size(), 'i')
        working_copy.lower_gaussian_elim(result)
        det = self.identity_element
        for i in range(self.rows):
            det = det * working_copy.data[i][i]
        return det

    def lup(self):
        """
        Function:       l, u, p = self.lup()
        Purpose:        Compute the LUP decomposition of self.
        Description:    This function returns three matrices
                        l, u, and p such that p * self = l * u
                        where l, u, and p have the following properties:

                        l is lower triangular with ones on the diagonal
                        u is upper triangular 
                        p is a permutation matrix.

                        The idea behind the algorithm is to first
                        do Gaussian elimination to obtain an upper
                        triangular matrix u and lower triangular matrix
                        r such that r * self = u, then by inverting r to
                        get l = r ^-1 we obtain self = r^-1 * u = l * u.
                        Note tha since r is lower triangular its
                        inverse must also be lower triangular.

                        Where does the p come in?  Well, with some
                        matrices our technique doesn't work due to
                        zeros appearing on the diagonal of r.  So we
                        apply some permutations to the original to
                        prevent this.
                        
        """
        upper = self.copy()
        result_inv = self.make_similar_matrix(self.size(), 'i')
        perm = self.make_similar_matrix((self.rows, self.rows), 'i')

        row_index, col_index = 0, 0
        last_row = self.rows - 1
        last_col = self.cols - 1
        while row_index < last_row and col_index < last_col:
            leader = upper.find_row_leader(row_index, col_index)
            if leader < 0:
                col_index = col_index+1
                continue
            if leader != row_index:
                upper.swap_rows(leader, row_index)
                result_inv.swap_rows(leader, row_index)
                perm.swap_rows(leader, row_index)
            row_index, col_index = upper.partial_lower_gauss_elim(
                row_index, col_index, result_inv)

        lower = self.make_similar_matrix((self.rows, self.rows), 'i')
        result_inv.lower_gaussian_elim(lower)
        result_inv.upper_inverse(lower)
        # possible optimization: due perm*lower explicitly without
        # relying on the * operator.
        return perm*lower, upper, perm

    def solve(self, b):
        """
        solve(self, b):

        b:   a list.

        Returns the values of x such that Ax = b.

        This is done using the LUP decomposition by 
        noting that Ax = b implies PAx = Pb implies LUx = Pb.
        First we solve for Ly = Pb and then we solve Ux = y.
        The following is an example of how to use solve:

>>> # Floating point example
>>> import genericmatrix
>>> A = genericmatrix.GenericMatrix(size=(2, 5), f_str=lambda x: '%.4f' % x)
>>> A.set_row(0, [0.0, 0.0, 0.160, 0.550, 0.280])
>>> A.set_row(1, [0.0, 0.0, 0.745, 0.610, 0.190])
>>> A
<matrix
 0.0000 0.0000 0.1600 0.5500 0.2800
 0.0000 0.0000 0.7450 0.6100 0.1900>
>>> b = [0.975, 0.350]
>>> x = A.solve(b)
>>> z = A.left_mul_column_vec(x)
>>> diff = reduce(lambda xx, yy: xx+yy, map(lambda aa, bb: abs(aa-bb), b, z))
>>> diff > 1e-6
0
>>> # Boolean example
>>> XOR = lambda x, y: x^y
>>> AND = lambda x, y: x&y
>>> DIV = lambda x, y: x
>>> m=GenericMatrix(size=(3, 6), zero_element=0, identity_element=1, add=XOR, mul=AND, sub=XOR, div=DIV)
>>> m.set_row(0, [1, 0, 0, 1, 0, 1])
>>> m.set_row(1, [0, 1, 1, 0, 1, 0])
>>> m.set_row(2, [0, 1, 0, 1, 1, 0])
>>> b = [0, 1, 1]
>>> x = m.solve(b)
>>> z = m.left_mul_column_vec(x)
>>> z
[0, 1, 1]

        """
        assert self.cols >= self.rows
        
        lower, upper, perm = self.lup()
        perm_b = perm.left_mul_column_vec(b)
        y = [0]*len(perm_b)
        for row in range(lower.rows):
            y[row] = perm_b[row]
            for i in range(row+1, lower.rows):
                perm_b[i] = lower.sub(perm_b[i], lower.mul(lower[i, row],
                                                           perm_b[row]))
        x = [0]*self.cols

        for cur_row in range(len(y)-1, -1, -1):
            col = upper.find_col_leader(cur_row, 0)
            assert col > -1
            x[col] = upper.div(y[cur_row], upper[cur_row, col])
            y[cur_row] = x[col]
            for i in range(0, cur_row):
                y[i] = upper.sub(y[i], upper.mul(upper[i, col], y[cur_row]))
        return x


def dot_product(mul, add, x, y):
    """
    Function:    dot_product(mul, add, x, y)
    Description: Return the dot product of lists x and y using mul and
                 add as the multiplication and addition operations.
    """
    assert len(x) == len(y), 'sizes do not match'
    return reduce(add, map(mul, x, y))


class GenericMatrixTester:
    def test(self, num_tests, size_list):
        """
        Function:       test(num_tests, size_list)

        Description:    For each test, run num_tests tests for square
                        matrices with the sizes in size_list.
        """

        for size in size_list:
            self.random_inverse_test(size, num_tests)
            self.random_lup_test(size, num_tests)
            self.random_solve_test(size, num_tests)
            self.random_det_test(size, num_tests)
            self.random_add_test(size, size + 1, num_tests)

    def make_random(self, s):
        import random 
        r = GenericMatrix(size=s, fill_mode=lambda x, y: random.random(),
                          equals_zero=lambda x: abs(x) < 1e-6)
        return r

    def mat_abs(self, m):
        r = -1
        rows, cols = m.size()
        for i in range(0, rows):
            for j in range(0, cols):
                if abs(m[i, j]) > r:
                    r = abs(m[i, j])
        return r

    def random_add_test(self, rows, cols, times):
        """Do some random tests to check addition and subtraction."""
        for i in range(times):
            first = self.make_random((rows, cols))
            second = self.make_random((rows, cols))
            m_sum = first + second
            assert self.mat_abs(m_sum - first - second) < 1e-6

    def random_inverse_test(self, s, n):
        ident = GenericMatrix(size=(s, s), fill_mode='i')
        for i in range(n):
            m = self.make_random((s, s))
            assert self.mat_abs(ident - m * m.inverse()) < 1e-6, (
                'offender = ' + repr(m))

    def random_lup_test(self, s, n):
        for i in range(n):
            m = self.make_random((s, s))
            lower, upper, perm = m.lup()
            assert self.mat_abs(perm * m - lower * upper) < 1e-6, \
                'offender = ' + repr(m)

    def random_solve_test(self, s, n):
        import random
        if s <= 1:
            return
        extra_equations = 3
        
        for i in range(n):
            m = self.make_random((s, s + extra_equations))
            for j in range(extra_equations):
                col_to_kill = random.randrange(s+extra_equations)
                for r in range(m.rows):
                    m[r, col_to_kill] = 0.0
            b = list(map(lambda x: random.random(), range(s)))
            x = m.solve(b)
            z = m.left_mul_column_vec(x)
            diff = reduce(lambda xx, yy: xx+yy,
                          map(lambda aa, bb: abs(aa-bb), b, z))
            assert diff < 1e-6, ('offenders: m = ' + repr(m) +
                                 '\nx = ' + repr(x) + '\nb = ' + repr(b) +
                                 '\ndiff = ' + repr(diff))

    def random_det_test(self, s, n):
        for i in range(n):
            m1 = self.make_random((s, s))
            m2 = self.make_random((s, s))
            prod = m1 * m2
            assert (abs(m1.determinant() * m2.determinant()
                        - prod.determinant())
                    < 1e-6), 'offenders = ' + repr(m1) + repr(m2)


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
The GenericMatrixTester class contains some simple
testing functions such as random_inverse_test, random_lup_test,
random_solve_test, and random_det_test which generate random floating
point values and test the appropriate routines.  The simplest way to
run these tests is via

>>> import genericmatrix
>>> t = genericmatrix.GenericMatrixTester()
>>> t.test(100, [1, 2, 3, 4, 5, 10])

# runs 100 tests each for sizes 1-5, 10
# note this may take a few minutes

If any problems occur, assertion errors are raised.  Otherwise
nothing is returned.  Note that you can also use the doctest
package to test all the python examples in the documentation
by typing 'python genericmatrix.py' or 'python -v genericmatrix.py' at the
command line.
"""


# The following code is used to make the doctest package
# check examples in docstrings when you enter
__test__ = {
    'testing_doc': testing_doc
}


def _test():
    import doctest
    import genericmatrix

    return doctest.testmod(genericmatrix)


if __name__ == "__main__":
    _test()
    print('Tests Finished.')
