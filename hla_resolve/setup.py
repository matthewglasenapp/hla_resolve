from setuptools import setup, find_packages

setup(
    name="hla_resolve",
    version="0.1.0",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "hla_resolve=hla_resolve.cli:main",
        ],
    },
)