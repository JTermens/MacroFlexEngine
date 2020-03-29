import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PPI",
    version='1.0',
    description='Protein-Protein Complex Builder',
    author='Miguel Luengo, Natalia Pattarone, Joan Termens',
    author_email='our emails - TBC',
    url='GitHub Project URL',
    py_modules=['README', 'launch', 'Alignment/superimposition'],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
