import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="MacroFlexEngine",
    version='1.0',
    description='Protein-Protein and Protein-DNA/RNA Complex Builder',
    author='Miguel Luengo, Natalia Pattarone, Joan Termens',
    author_email='our emails - TBC',
    url='GitHub Project URL',
    py_modules=['README', 'MFEngine-launch', 'Alignment/superimposition'],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "biopython >= 1.73",
        "numpy >= 1.16.2",
    ],
    python_requires='>=3.6',
)