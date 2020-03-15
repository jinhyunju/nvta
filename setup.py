from setuptools import setup, find_packages

setup(name="nvta",
      version='0.1',
      description="Mapping transcript coordinates to the reference genome based on their CIGAR string",
      author="Jin Hyun Ju",
      author_email="jinhyun.ju@gmail.com",
      packages= ['nvta'],
      install_requires=["intervaltree", "pytest"]
      )
