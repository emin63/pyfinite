"""Some basic and quick unit tests.
"""

import unittest
import operator
import types

from pyfinite import ffield
from pyfinite import genericmatrix


class MatrixTester(unittest.TestCase):
    """Class to test various matrix operations.
    """

    def do_simple_determinant_test(self, my_field, element, answer):        
        """Do a simple test of determinant calculation.

        :param my_field:   Instance of FField or at leat something with the
                           Add, Subtract, Multiply, and Divide methods.

        :param element:    Field element to use in calculations.

        :param answer:     Expected answer.

        ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

        PURPOSE:  Construct a simple square, 2x2 matrix and check
                  that the determinant is correct.

        """
        matrix = genericmatrix.GenericMatrix(
            (2, 2), zeroElement=0, identityElement=1,
            add=my_field.Add, sub=my_field.Subtract, mul=my_field.Multiply,
            div=my_field.Divide)
        matrix.SetRow(0, [element, 0])
        matrix.SetRow(1, [0, element])
        self.assertEqual(matrix.Determinant(), answer)

    def test_determinants(self):
        """Test that the determinant is calculated correctly.
        """

        # First we check F16 where we use the shortcut of passing in
        # integers for the field elements.
        my_field = ffield.FField(4)  #  F16
        self.do_simple_determinant_test(my_field, 15, 10)

        # Next we try the same thing where we explicitly pass in the
        # field elements.
        self.do_simple_determinant_test(my_field, ffield.FElement(
            my_field, 15), ffield.FElement(my_field, 10))

        # Next we try the same thing with reals
        my_field = types.SimpleNamespace(
            Add=operator.add, Subtract=operator.sub, Multiply=operator.mul,
            Divide=operator.truediv)  # regular field of reals
        self.do_simple_determinant_test(my_field, 15, 225)
            
