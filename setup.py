# read the contents of your README file
from pathlib import Path

from setuptools import setup

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="scubas",
    version="0.1.5",
    packages=["scubas"],
    package_dir={"scubas": "scubas"},
    package_data={"scubas": []},
    author="Shibaji Chakraborty",
    author_email="shibaji7@vt.edu",
    maintainer="Shibaji Chakraborty",
    maintainer_email="shibaji7@vt.edu",
    license="MIT License",
    description="SCUBAS: Submarine Cables Upset by Auroral Streams.",
    # long_description="Model compute electrical surges in submarine cables induced by geomagnetic activities",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=[
        "pandas",
        "matplotlib>=3.3.4",
        "pyproj",
        "loguru",
        "scipy",
        "SciencePlots",
    ],
    keywords=[
        "Python>=3.6",
        "submarine cable",
        "electrical surge",
        "geomagnetic induction",
        "geomagnetic disturbance",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Education",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Operating System :: OS Independent",
    ],
    project_urls={
        "Documentations": "https://scubas.readthedocs.io/en/latest/",
    },
    url="https://github.com/shibaji7/SCUBAS",
)
