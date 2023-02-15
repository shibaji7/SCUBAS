from setuptools import setup

setup(
    name="scubas",
    version="0.1.0",
    packages=["scubas"],
    package_dir={"scubas": "scubas"},
    package_data={"scubas": []},
    author="Shibaji Chakraborty",
    author_email="shibaji7@vt.edu",
    maintainer="Shibaji Chakraborty",
    maintainer_email="shibaji7@vt.edu",
    license="MIT License",
    description="SCUBAS: Submarine Cables Upset by Auroral Streams.",
    long_description="Model compute electrical surges in submarine cables induced by geomagnetic activities",
    install_requires=[
        "pandas",
        "matplotlib>=3.2",
    ],
    keywords=[
        "python",
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
)
