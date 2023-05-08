from setuptools import setup
from setuptools import Command
from setuptools import Extension
from Cython.Build import cythonize

try:
    import sage.env
except ImportError:
    raise ValueError("this package requires SageMath")

class TestCommand(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import subprocess
        if subprocess.call(['sage', '-tp', '--force-lib', 'src/']):
            raise SystemExit("Doctest failures")

def do_cythonize():
    return cythonize(
            [Extension(
                "*",
                ["src/operator_gb/rational_linear_algebra.pyx"],
                extra_compile_args=['-std=c11'],
            )],
            aliases = sage.env.cython_aliases(),
        )

try:
    from sage.misc.package_dir import cython_namespace_package_support
    with cython_namespace_package_support():
        extensions = do_cythonize()
except ImportError:
    extensions = do_cythonize()

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name = "operator_gb",
    version = "1.0",
    author = "Clemens Hofstadler",
    author_email = "clemens@hofstadler@mathematik.uni-kassel.de",
    license = "GPLv2",
    packages = [
        "operator_gb",
    ],
    install_requires = required,
    package_dir = {'': 'src/'},
    ext_modules = extensions,
    include_dirs = sage.env.sage_include_directories() + ["."],
    cmdclass = {'test': TestCommand},
    zip_safe=False,
)