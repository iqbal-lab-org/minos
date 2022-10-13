from setuptools import setup, find_packages
import shutil


if shutil.which("gramtools") is None:
    raise RuntimeError("Please install gramtools. Cannot continue")

with open("requirements.txt") as f:
    install_requires = [x.rstrip() for x in f]

setup(
    name="bio-minos",
    version="0.12.5",
    description="Variant call adjudication",
    packages=find_packages(),
    package_data={"minos": ["test_data_files/*"]},
    author="Martin Hunt",
    author_email="mhunt@ebi.ac.uk",
    url="https://github.com/iqbal-lab-org/minos",
    test_suite="nose.collector",
    tests_require=["pytest"],
    entry_points={"console_scripts": ["minos = minos.__main__:main"]},
    install_requires=install_requires,
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: MIT License",
    ],
)
