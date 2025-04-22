import glob

import setuptools
import yaml

# Version
version = ""
name = ""
fversion = glob.glob("**/_version.py", recursive=True)[0]
with open(fversion) as fid:
    lines = fid.read().splitlines()
    name = lines[0].split("=")[-1].strip().replace('"', "")
    version = lines[1].split("=")[-1].strip().replace('"', "")

# App name - dependencies
env = {}
"""
with open("recipes/workflow.yaml") as fid:
    env = yaml.safe_load(fid)
"""
install_requires = []
"""
for package in env["dependencies"]:
    if isinstance(package, dict):
        package = package.get("pip", [])
        install_requires += package
    else:
        install_requires.append(package)
"""
description = "R-group expansion of ChEBI molecules"

setuptools.setup(
    name=name,
    version=version,
    author=["guillaume-gricourt"],
    author_email=["guipagui@gmail.com"],
    description=description,
    long_description_content_type="text/markdown",
    url="https://github.com/brsynth/chebi-rgroup",
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    include_package_data=True,
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    install_requires=install_requires,
)
