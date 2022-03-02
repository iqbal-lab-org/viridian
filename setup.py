import glob
from setuptools import setup, find_packages


with open("requirements.txt") as f:
    install_requires = [x.rstrip() for x in f]

setup(
    name="viridian_workflow",
    version="0.3.7",
    description="FIXME",
    packages=find_packages(),
    package_data={"viridian_workflow": ["amplicon_scheme_data/*"]},
    author="Jeff Knaggs,Martin Hunt",
    author_email="FIXME",
    url="https://github.com/iqbal-lab-org/viridian_workflow",
    test_suite="nose.collector",
    tests_require=["pytest"],
    entry_points={
        "console_scripts": ["viridian_workflow = viridian_workflow.__main__:main"]
    },
    install_requires=install_requires,
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: MIT License",
    ],
)
