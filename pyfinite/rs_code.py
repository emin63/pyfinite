
# Copyright Emin Martinian 2002.  See below for license terms.
# Version Control Info: $Id: rs_code.py,v 1.4 2008-01-05 22:08:45 emin Exp $

"""
This package implements the RSCode class designed to do
Reed-Solomon encoding and (erasure) decoding.  The following
docstrings provide detailed information on various topics.

  RSCode.__doc__   Describes the RSCode class and how to use it.

  license_doc      Describes the license and lack of warranty.

"""

import ffield
import genericmatrix
import math
import doctest


class RSCode:
    """
    The RSCode class implements a Reed-Solomon code
    (currently, only erasure decoding not error decoding is
    implemented).  The relevant methods are:

    __init__
    encode
    decode_immediate
    decode
    prepare_decoder
    random_test

    A breif example of how to use the code follows:

>>> import rs_code

# Create a coder for an n, k = 16, 8 code and test
# decoding for a simple erasure pattern.

>>> C = rs_code.RSCode(16,8)
>>> inVec = list(range(8))
>>> codedVec = C.encode(inVec)
>>> receivedVec = list(codedVec)

# now erase some entries in the encoded vector by setting them to None
>>> receivedVec[3] = receivedVec[9] = receivedVec[12] = None
>>> receivedVec
[0, 1, 2, None, 4, 5, 6, 7, 8, None, 10, 11, None, 13, 14, 15]
>>> decVec = C.decode_immediate(receivedVec)
>>> decVec
[0, 1, 2, 3, 4, 5, 6, 7]

# Now try the random testing method for more complete coverage.
# Note this will take a while.
>>> for k in range(1, 8):
...     for p in range(1, 12):
...         c = rs_code.RSCode(k+p, k)
...         c.random_test(25)
>>> for k in range(1, 8):
...     for p in range(1, 12):
...         c = rs_code.RSCode(k+p, k, systematic=0)
...         c.random_test(25)
"""

    def __init__(self, n, k, log_2_field_size=-1, systematic=1,
                 should_ese_lut=-1):
        """
        Function:  __init__(n, k, log_2_field_size, systematic, should_ese_lut)
        Purpose:   Create a Reed-Solomon coder for an (n, k) code.
        Notes:     The last parameters, log_2_field_size, systematic
                   and should_ese_lut are optional.

                   The log_2_field_size parameter
                   represents the base 2 logarithm of the field size.
                   If it is omitted, the field GF(2^p) is used where
                   p is the smallest integer where 2^p >= n.

                   If systematic is true then a systematic encoder
                   is created (i.e. one where the first k symbols
                   of the encoded result always match the data).

                   If should_ese_lut == 1 then a lookup table is used for
                   computing finite field multiplies and divides.
                   If should_ese_lut == 0 then no lookup table is used.
                   If should_ese_lut == -1 (the default), then the code
                   decides when a lookup table should be used.
        """
        if log_2_field_size < 0:
            log_2_field_size = int(math.ceil(math.log(n) / math.log(2)))
        self.field = ffield.FField(log_2_field_size, use_lut=should_ese_lut)
        self.n = n
        self.k = k
        self.fieldSize = 1 << log_2_field_size
        self.encoder_matrix = None
        self.decoder_matrix = None
        self.create_encoder_matrix()
        if systematic:
            self.encoder_matrix.transpose()
            self.encoder_matrix.lower_gaussian_elim()
            self.encoder_matrix.upper_inverse()
            self.encoder_matrix.transpose()

    def __repr__(self):
        rep = ('<RSCode (n, k) = (' + repr(self.n) + ', ' + repr(self.k) + ')'
               + '  over GF(2^' + repr(self.field.p) + ')\n' +
               repr(self.encoder_matrix) + '\n' + '>')
        return rep

    def create_encoder_matrix(self):
        self.encoder_matrix = genericmatrix.GenericMatrix(
            (self.n, self.k), 0, 1, self.field.add, self.field.subtract,
            self.field.multiply, self.field.divide)
        self.encoder_matrix[0, 0] = 1
        for i in range(0, self.n):
            term = 1
            for j in range(0, self.k):
                self.encoder_matrix[i, j] = term
                term = self.field.multiply(term, i)

    def encode(self, data):
        """
        Function:       encode(data)
        Purpose:        encode a list of length k into length n.
        """
        assert len(data) == self.k, 'encode: input data must be size k list.'

        return self.encoder_matrix.left_mul_column_vec(data)

    def prepare_decoder(self, un_erased_locations):
        """
        Function:       prepare_decoder(un_erased_locations)
        Description:    The input un_erased_locations is a list of the first
                        self.k elements of the codeword which were
                        NOT erased.  For example, if the 0th, 5th,
                        and 7th symbols of a (16, 5) code were erased,
                        then prepare_decoder([1, 2, 3, 4, 6]) would
                        properly prepare for decoding.
        """
        if len(un_erased_locations) != self.k:
            raise ValueError('input must be exactly length k')

        limited_encoder = genericmatrix.GenericMatrix(
            (self.k, self.k), 0, 1, self.field.add, self.field.subtract,
            self.field.multiply, self.field.divide)
        for i in range(0, self.k):
            limited_encoder.set_row(
                i, self.encoder_matrix.get_row(un_erased_locations[i]))
        self.decoder_matrix = limited_encoder.inverse()

    def decode(self, un_erased_terms):
        """
        Function:       decode(un_erased_terms)
        Purpose:        Use the
        Description:
        """
        return self.decoder_matrix.left_mul_column_vec(un_erased_terms)

    def decode_immediate(self, data):
        """
        Function:       decode_immediate(data)
        Description:    Takes as input a data vector of length self.n
                        where erased symbols are set to None and
                        returns the decoded result provided that
                        at least self.k symbols are not None.

                        For example, for an n, k = 6, 4 code, a
                        decodable input vector would be
                        [2, 0, None, 1, 2, None].
        """

        if len(data) != self.n:
            raise ValueError('input must be a length n list')

        un_erased_locations = []
        un_erased_terms = []
        for i in range(self.n):
            if data[i] is not None:
                un_erased_locations.append(i)
                un_erased_terms.append(data[i])
        self.prepare_decoder(un_erased_locations[0:self.k])
        return self.decode(un_erased_terms[0:self.k])

    def random_test(self, num_tests):
        import random

        max_erasures = self.n-self.k
        for i in range(num_tests):
            in_vec = list(range(self.k))
            for j in range(self.k):
                in_vec[j] = random.randint(0, (1 << self.field.p) - 1)
            coded_vec = self.encode(in_vec)
            num_erasures = random.randint(0, max_erasures)
            for j in range(num_erasures):
                j = random.randint(0, self.n-1)
                while coded_vec[j] is None:
                    j = random.randint(0, self.n-1)
                coded_vec[j] = None
            dec_vec = self.decode_immediate(coded_vec)
            assert dec_vec == in_vec, ('in_vec = ' + repr(in_vec)
                                       + '\ncoded_vec = ' + repr(coded_vec)
                                       + '\ndec_vec = ' + repr(dec_vec))


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


# The following code is used to make the doctest package
# check examples in docstrings.
def _test():
    return doctest.testmod()


if __name__ == "__main__":
    _test()
    print('Tests Finished')
