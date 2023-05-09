# operator_gb

## Description

SageMath package for Gröbner basis computations in the free algebra, with dedicated methods for automatising the proofs of operator statements.
In particular, the package provides several methods for searching for elemnts of certain forms in noncommutative polynomial ideals.

## License

Distributed under the terms of the GNU General Public License (GPL, see the
LICENSE file), either version 2 or (at your option) any later version

- http://www.gnu.org/licenses/

## Requirements

- SageMath 9.1 or later is recommended

## Dependencies

- The Python library `pyahocorasick` (https://pyahocorasick.readthedocs.io/en/latest/)

## Installation

### With SageMath built from source or binaries from sagemath.org

**Note**: This way of installing the package also automatically installs the `pyahocorasick` library.
Thus, except executing the command below, no additional work is required.

To download and install the latest version on a system where SageMath
was built from source or installed from official packages, run

    sage -pip install [--user] git+https://github.com/ClemensHofstadler/operator_gb.git
The optional `--user` flag causes the package to be installed in your `.sage`
directory instead of the SageMath installation tree.

Alternatively, run (square brackets indicate optional flags)

    sage -pip install [--user] [--editable] .

from the root of a local git checkout. The `--editable` flag causes the
"installed" version to point to your local checkout, making it easier,
if you edit the code, to run the modified version. See the pip documentation
for more installation options.

Microsoft Windows users should run the above commands in a "SageMath shell", see

- https://wiki.sagemath.org/SageWindows

Apple macOS users may need additional steps before they are able to add external
packages to their Sage installations. See

- https://github.com/3-manifolds/fix_mac_sage/releases
- https://ask.sagemath.org/question/51130

for more information.

### Using the package without installation

**Note**: This way of using the package requires the library `pyahocorasick` to be 
already installed (or at least visible to Sage).
 
To use operator_gb directly from a git checkout (without installation), run

    sage -python setup.py build_ext --inplace

from the checkout, and add the `src/` directory to your Python `sys.path`.

The package contains compiled (Cython) modules which are automatically built as
part of the installation procedure. Installation will fail if they cannot be
built.


## Contact

Clemens Hofstadler (<clemens.hofstadler@mathematik.uni-kassel.de>)

