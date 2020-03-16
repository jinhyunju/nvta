import sys

from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand


class PyTest(TestCommand):

    # Adopted from
    # https://pytest.readthedocs.io/en/2.7.3/goodpractises.html
    user_options = []

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        # import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)


setup(name="nvta",
      version='0.1',
      description=("Mapping transcript coordinates to the "
                   "reference genome based on their CIGAR string"),
      author="Jin Hyun Ju",
      author_email="jinhyun.ju@gmail.com",
      packages=['nvta'],
      install_requires=["intervaltree"],
      test_suite="pytest",
      tests_require=["pytest"],
      cmdclass={'test': PyTest}
      )
