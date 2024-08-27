from setuptools import setup, find_packages

setup(
    name="mappability-diff",
    version="0.1",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
)