
# Copyright Emin Martinian 2002.  See below for license terms.
# Version Control Info: $Id: file_ecc.py,v 1.5 2008-01-05 22:08:44 emin Exp $

__doc__ = """
This package implements an erasure correction code for files.
Specifically it lets you take a file F and break it into N
pieces (which are named F.p_0, F.p_1, ..., F.p_N-1) such that
F can be recovered from any K pieces.  Since the size of each
piece is F/K (plus some small header information).

How is this better than simply repeating copies of a file?

Firstly, this package lets you get finer grained
redundancy control since producing a duplicate copy of a file
requires at least 100% redundancy while this package lets you
expand the redundancy by n/k (e.g. if n=11, k=10 only 10%
redundancy is added).

Secondly, using a Reed-Solomon code as is done in this package,
allows better loss resistance.  For example, assume you just
divided a file F into 4 pieces, F.1, F.2, ..., F.4, duplicated
each piece and placed them each on a different disk.  If the
two disks containing a copy of F.1 go down then you can no longer
recover F.

With the Reed-Solomon code used in this package, if you use n=8, k=4
you divide F into 8 pieces such that as long as at least 4 pieces are
available recovery can occur.  Thus if you placed each piece on a
separate disk, you could recover data as if any combination of 4 or
less disks fail.

The docstrings for the functions encode_file and decode_files
provide detailed information on usage and the docstring
license_doc describes the license and lack of warranty.

The following is an example of how to use this file:

>>> import file_ecc
>>> test_file = '/bin/ls'      # A reasonable size file for testing.
>>> prefix = '/tmp/ls_backup'  # Prefix for shares of file.
>>> # break into N=15 pieces
>>> names = file_ecc.encode_file(test_file, prefix, 15, 11)

# Imagine that only pieces 0, 1, 5, 4, 13, 8, 9, 10, 11, 12 & 14 are available.
>>> dec_list = [prefix + '.p_' + repr(x) for x in [0,1,5,4,13,8,9,10,11,12,14]]

>>> decoded_file = '/tmp/ls.r'  # Choose where we want reconstruction to go.
>>> file_ecc.decode_files(dec_list, decoded_file)
>>> fd1 = open(test_file,'rb')
>>> fd2 = open(decoded_file,'rb')
>>> fd1.read() == fd2.read()
1
"""


import os
import struct
from array import array

from rs_code import RSCode

headerSep = '|'


def get_file_size(filename):
    return os.stat(filename)[6]


def make_header(filename, n, k, size):
    return headerSep.join([
        'RS_PARITY_PIECE_HEADER', 'FILE', filename,
        'n', repr(n), 'k', repr(k), 'size', repr(size), 'piece']) + headerSep


def parse_header(header):
    if isinstance(header, bytes):
        header = header.decode('UTF8')
    return header.split(headerSep)


def read_encode_and_write_block(read_size, in_fd, out_fd, code):
    buffer = array('B')
    buffer.fromfile(in_fd, read_size)
    for i in range(read_size, code.k):
        buffer.append(0)
    code_vec = code.encode(buffer)
    for j in range(code.n):
        out_fd[j].write(struct.pack('B', code_vec[j]))


def encode_file(filename, prefix, n, k):
    """
    Function:     encode_file(filename, prefix, n, k)
    Description:  Encodes the file named by filename into n pieces named
                  prefix.p_0, prefix.p_1, ..., prefix.p_n-1.  At least
                  k of these pieces are needed for recovering filename.
                  Each piece is roughly the size of filename / k (there
                  is a very small overhead due to some header information).

                  Returns a list containing names of files for the pieces.

                  Note n and k must satisfy 0 < k < n < 257.
                  Use the decode_files function for decoding.
    """
    file_list = []
    if n > 256 or k >= n or k <= 0:
        raise Exception('Invalid (n, k), need 0 < k < n < 257.')
    in_fd = open(filename, 'rb')
    in_size = get_file_size(filename)
    header = make_header(filename, n, k, in_size)
    code = RSCode(n, k, 8, shouldUseLUT=-(k != 1))
    out_fd = list(range(n))
    for i in range(n):
        out_file_name = prefix + '.p_' + repr(i)
        file_list.append(out_file_name)
        out_fd[i] = open(out_file_name, 'wb')
        out_fd[i].write((header + repr(i) + '\n').encode('UTF8'))

    if k == 1:  # just doing repetition coding
        s = in_fd.read(1024)
        while s:
            _ = [x.write(s) for x in out_fd]
            s = in_fd.read(256)
    else:  # do the full blown RS encoding
        for i in range(0, (in_size//k)*k, k):
            read_encode_and_write_block(k, in_fd, out_fd, code)
        
        if (in_size % k) > 0:
            read_encode_and_write_block(in_size % k, in_fd, out_fd, code)

    return file_list


def extract_piece_nums(filenames, headers):
    li = list(range(len(filenames)))
    piece_nums = list(range(len(filenames)))
    for i in range(len(filenames)):
        li[i] = parse_header(headers[i])
    for i in range(len(filenames)):
        if (li[i][0] != 'RS_PARITY_PIECE_HEADER' or
                li[i][2] != li[0][2] or li[i][4] != li[0][4] or
                li[i][6] != li[0][6] or li[i][8] != li[0][8]):
            raise Exception(
                'File ' + repr(filenames[i]) + ' has incorrect header.')
        piece_nums[i] = int(li[i][10])
    n, k, size = int(li[0][4]), int(li[0][6]), int(li[0][8])
    if len(piece_nums) < k:
        raise Exception(('Not enough parity for decoding; needed ' +
                         repr(li[0][6])+' got '+repr(len(filenames))+'.'))
    return n, k, size, piece_nums


def read_decode_and_write_block(write_size, in_fds, out_fd, code):
    buffer = array('B')
    for j in range(code.k):
        buffer.fromfile(in_fds[j], 1)
    result = code.decode(buffer.tolist())
    for j in range(write_size):
        out_fd.write(struct.pack('B', result[j]))


def decode_files(filenames, out_name):
    """
    Function:     decode_files(filenames, out_name)
    Description:  Takes pieces of a file created using encode_files and
                  recovers the original file placing it in outName.
                  The argument filenames must be a list of at least k
                  file names generated using encode_files.
    """
    in_fds = list(range(len(filenames)))
    headers = list(range(len(filenames)))
    for i in range(len(filenames)):
        in_fds[i] = open(filenames[i], 'rb')
        headers[i] = in_fds[i].readline()
    n, k, in_size, piece_nums = extract_piece_nums(filenames, headers)
    out_fd = open(out_name, 'wb')
    code = RSCode(n, k, 8)
    dec_list = piece_nums[0:k]
    code.prepare_decoder(dec_list)
    for i in range(0, (in_size//k)*k, k):
        read_decode_and_write_block(k, in_fds, out_fd, code)
    if (in_size % k) > 0:
        read_decode_and_write_block(in_size % k, in_fds, out_fd, code)


license_doc = """
  This code was originally written by Emin Martinian (emin@alum.mit.edu).
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
    import doctest
    import file_ecc

    return doctest.testmod(file_ecc)


if __name__ == "__main__":
    _test()
    print('Tests finished')
