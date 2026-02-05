import os
import subprocess
import shutil

from setuptools import setup
from setuptools import Command
from setuptools import Extension
from setuptools.command.build import build
from setuptools.command.install import install
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

class BuildWithCadical(build):
    def run(self):
        self._build_cadical()
        super().run()

    def _build_cadical(self):
        cadical_dir = os.path.join(
            os.path.dirname(__file__),
            "includes",
            "cadical",
        )

        # Run configure and make
        print("Configuring CaDiCaL...")
        subprocess.check_call(["./configure"], cwd=cadical_dir)
        print("Building CaDiCaL...")
        subprocess.check_call(["make"], cwd=cadical_dir)
        
class InstallWithCadical(install):
    """Copy the built CaDiCaL binary into the installed Python package."""
    def run(self):
        super().run()
        self._install_cadical()

    def _install_cadical(self):
        cadical_bin = os.path.join(
            os.path.dirname(__file__),
            "includes",
            "cadical",
            "build",
            "cadical"
        )

        if not os.path.isfile(cadical_bin):
            raise RuntimeError("CaDiCaL binary not found; build might have failed")

        # Copy into the installed package directory
        package_dir = os.path.join(self.install_lib, "operator_gb")
        os.makedirs(package_dir, exist_ok=True)
        dst = os.path.join(package_dir, "cadical")
        shutil.copy2(cadical_bin, dst)
        os.chmod(dst, 0o755)

def do_cythonize():
    return cythonize(
            [Extension(
                "*",
                ["src/operator_gb/*.pyx"],
            )],
            aliases = sage.env.cython_aliases(required_modules=(), optional_modules=()),
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
    version = "1.1",
    author = "Clemens Hofstadler",
    author_email = "clemens@hofstadler@jku.at",
    license = "GPLv2",
    packages = [
        "operator_gb",
    ],
    install_requires = required,
    package_dir = {'': 'src/'},
    ext_modules = extensions,
    include_dirs=sage.env.sage_include_directories() + ["."],
  cmdclass={
    'test': TestCommand,
    'build': BuildWithCadical,
    'install': InstallWithCadical,
    },
    zip_safe=False,
)
