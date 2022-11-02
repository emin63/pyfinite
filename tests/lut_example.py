"""Simple examples of doctests to illustrate lookup table (LUT).

You can run theses tests interactively in a python interpreter
or you can do something like

    python3 <PATH_TO_THIS_FILE>

on the command line to test everything.

To make sure we start from a clean environment, you will first
want to make sure you delete any existing lookup tables.

>>> import os
>>> from pyfinite import ffield
>>> lutName = ffield.FField.make_lut_path(8)
>>> if os.path.exists(lutName):
...     os.remove(lutName)
...


If you are trying to debug something, you may want to explicitly
disable the lookup table by passing `useLUT=0` when creating a field.
That will prevent you from using an old lookup table built from
a different generator.  The following shows an example where we
disable the lookup table:

>>> from pyfinite import ffield
>>> f_0x11B = ffield.FField(8, gen=0x11B, useLUT=0)
>>> res = f_0x11B.Multiply(0xFF, 2) 
>>> print(res)
229

Usually, lookup tables are built for small fields.

>>> from pyfinite import ffield
>>> f = ffield.FField(8)
>>> print(f.generator)
285
>>> print(f.Multiply(0xFF, 2))
227

If you now create a field using a lookup table, you may get
unexpected results:

>>> from pyfinite import ffield
>>> try:  # the following will raise an exception:
...     f_0x11B_with_LUT = ffield.FField(8, gen=0x11B)
... except ValueError as problem:
...     print(problem)
...
Refusing to use lookup table in file ffield.lut.8.
That file is for a different or unknown generator.
Please remove that file or pass useLUT=0 to init.


Another altrenative is to explicitly remove the lookup table
if it exists to prevent this issue:

>>> import os
>>> from pyfinite import ffield
>>> lutName = ffield.FField.make_lut_path(8)
>>> if os.path.exists(lutName):
...     os.remove(lutName)
...
>>> f_0x11B_with_LUT = ffield.FField(8, gen=0x11B)
>>> res = f_0x11B_with_LUT.Multiply(0xFF, 2) 
>>> print(res)
229


"""

import doctest


if __name__ == '__main__':
    doctest.testmod()
    print('Finished Tests')
