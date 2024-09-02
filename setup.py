from pathlib import Path
from setuptools import setup, find_packages

author = "Jonas Simon Fleck, Zhisong He, Leander Dony"
author_email = "jonas.simon.fleck@gmail.com, zhisong.he@bsse.ethz.ch, leander.dony@helmholtz-munich.de"
description = "Human Neural Organoid Cell Atlas Toolbox"

long_description = Path("README.md").read_text("utf-8")
requirements = [
    l.strip() for l in Path("requirements.txt").read_text("utf-8").splitlines()
]

setup(
    name="hnoca",
    version="0.1.1",
    author=author,
    author_email=author_email,
    description=description,
    license="MIT",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    install_requires=requirements,
    python_requires=">=3.9",
)
