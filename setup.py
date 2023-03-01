"""A setuptools based setup module.

See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_namespace_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name="SCaFoSpy",
    version="0.0.1",
    description="..Description coming" ,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/iTaxoTools/SCaFoSpy/",
    author="David Leisse",
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3 :: Only",
    ],
    package_dir={"": "src"},
    packages=find_namespace_packages(
        include=("itaxotools*",),
        where="src",
    ),
    python_requires=">=3.9, <4",
    install_requires=[
        "Dendropy",
    ],
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "scafospy= itaxotools.SCaFoSpy.SCaFoS:run"
        ],
    #     "pyinstaller40": [
    #         "hook-dirs = itaxotools.__pyinstaller:get_hook_dirs",
    #         "tests = itaxotools.__pyinstaller:get_pyinstaller_tests",
    #     ],
    }
)
