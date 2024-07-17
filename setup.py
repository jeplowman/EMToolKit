from setuptools import setup, find_packages

setup(
    name="EMToolKit",
    version="0.1.0",
    author="Joseph Plowman",
    author_email="joseph.plowman@swri.org",
    description="A Standardized Framework for Computing and Visualizing Differential Emission Measures",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url="https://github.com/jeplowman/EMToolKit",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: BSD 3-Clause License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        "numpy>=1.19.2",
        "pandas>=1.1.3",
        "scipy>=1.5.2",
    ],
)