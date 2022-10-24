from setuptools import setup

setup(name='siRNAenricher',
      version='0.1.1',
      description='siRNAenricher',
      author='Stephen Fletcher',
      author_email='s.fletcher@uq.edu.au',
      license='MIT',
      packages=['siRNAenricher'],
      classifiers=[
    # How mature is this project? Common values are
    #   3 - Alpha
    #   4 - Beta
    #   5 - Production/Stable
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Science/Research',
    # Indicate who your project is intended for

    'Topic :: Scientific/Engineering :: Bio-Informatics',

    # Pick your license as you wish (should match "license" above)
     'License :: OSI Approved :: MIT License',

    # Specify the Python versions you support here. In particular, ensure
    # that you indicate whether you support Python 2, Python 3 or both.
    'Programming Language :: Python :: 3.11'],
      install_requires=['numpy'],
      zip_safe=False)