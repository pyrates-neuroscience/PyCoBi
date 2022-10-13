from setuptools import setup, find_packages

with open("requirements.txt", "r", encoding="utf8") as fh:
    REQUIREMENTS = fh.read()
INSTALL_REQUIREMENTS = REQUIREMENTS

CLASSIFIERS = ["Programming Language :: Python :: 3",
               "Programming Language :: Python :: 3.6",
               "Programming Language :: Python :: 3.7",
               "Programming Language :: Python :: 3.8",
               "Programming Language :: Python :: 3.9",
               "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
               "Operating System :: OS Independent",
               "Development Status :: 3 - Alpha",
               "Intended Audience :: Developers",
               "Intended Audience :: Science/Research",
               "Topic :: Scientific/Engineering",
               ]

EXTRAS = {"dev": ["pytest", "bump2version"]}

with open("VERSION", "r", encoding="utf8") as fh:
    VERSION = fh.read().strip()

with open("README.md", "r", encoding="utf8") as fh:
    DESCRIPTION = fh.read()

setup(name='pycobi',
      version=VERSION,
      description='Python tool for parameter continuation and bifurcation analysis',
      long_description=DESCRIPTION,
      long_description_content_type='text/markdown',
      author="Richard Gast",
      author_email='richard.gast@northwestern.edu',
      license='GPL v3',
      packages=find_packages(),
      zip_safe=False,
      python_requires='>=3.6',
      install_requires=INSTALL_REQUIREMENTS,
      extras_require=EXTRAS,
      classifiers=CLASSIFIERS,
      include_package_data=True  # include additional non-python files specified in MANIFEST.in
      )
