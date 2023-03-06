from setuptools import setup, find_packages


with open("requirements.txt") as f:
    install_requires = [x.rstrip() for x in f]

setup(
    name="viridian",
    version="0.4.0",
    description="Consensus builder from amplicon sequenced virus reads",
    packages=find_packages(),
    package_data={"viridian": ["amplicon_scheme_data/*"]},
    author="Jeff Knaggs,Martin Hunt",
    author_email="FIXME",
    url="https://github.com/iqbal-lab-org/viridian",
    tests_require=["pytest"],
    entry_points={
        "console_scripts": [
            "viridian = viridian.__main__:main",
            "viridian_workflow = viridian.__main__:main",
        ]
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
