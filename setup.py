import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pycctyper",
    version="1.0.0",
    author="Sam Nooij",
    author_email="s.nooij@uu.nl",
    description="PyCCTyper: Accelerated detection and subtyping of CRISPR-Cas",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/samnooij/PyCCTyper",
    download_url="https://github.com/samnooij/PyCCTyper/archive/v1.0.0.tar.gz",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Development Status :: 4 - Beta",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy >= 1.17.5, < 2",
        "pandas >= 1.3, < 3",
        "scipy >= 1.4.1, < 2",
        "biopython >= 1.76, < 2",
        "multiprocess >= 0.70.9, < 1",
        "scikit-learn >= 0.22.0, < 2",
        "xgboost >= 1.4, < 2",
        "tqdm >= 4, < 5",
        "drawSvg >= 1.8.0, < 2",
        "setuptools",
        "pyrodigal >= 3.7.1, < 4",
        "pyhmmer >= 0.12.1, < 1",
    ],
    scripts=["bin/pycctyper", "bin/repeatType", "bin/repeatTrain"],
)
