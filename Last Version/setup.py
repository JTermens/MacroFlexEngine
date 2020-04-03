import setuptools

with open("../README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="MacroFlexEngine",
    version='1.0',
    description='MacroFlexEngine is a Protein-Protein and Protein-DNA/RNA Complex Builder',
    author='Miguel Luengo, Natalia Pattarone, Joan Termens',
    author_email='miguel.luengo01@estudiant.upf.edu, natalia.pattarone01@estudiant.upf.edu and joan.termensi01@estudiant.upf.edu',
    url='GitHub Project URL',
    py_modules=['README', 'MFEngine-launch'],
    packages=['MacroFlexEngine','MacroFlexEngine/lib','MacroFlexEngine/Modeller'],
    scripts = ['MFEngine-launch'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "biopython >= 1.73",
        "numpy >= 1.16.2",
        "modeller",
        "setuptools",
        "regex",
        "matplotlib"
    ],
    python_requires='>=3.6',
)