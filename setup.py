from setuptools import find_packages, setup
from src.version import __version__

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="NucFreq",
    version=__version__,
    description="",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    author="Andrew Harris",
    author_email="ajharris@cvm.tamu.edu",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.10",
        "Operating System :: OS Independent",
    ],
    scripts=["src/nucfreq/__main__.py"],
    entry_points={"console_scripts": ["nucfreq=nucfreq.__main__:main"]},
    python_requires=">=3.10",
)
