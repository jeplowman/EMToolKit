# Installation Guide for EMToolKit

## PyPI Installation

EMToolKit is available on the [Python Package Index (PyPI)](https://pypi.org/project/EMToolKit/). It is strongly recommended to use a Python virtual environment:

### Installation Steps

Install and test in the shell:
``` sh
    python -m venv .venv
    source .venv/bin/activate
    pip install emtoolkit
    emtest
```

You can also test within Python:
``` python
    import EMToolKit.examples as ex
    ex.example_dir()
    ex.example_run()
```

## GitHub Installation
For the latest development version, clone the repository from GitHub and install it:
```
git clone https://github.com/jeplowman/EMToolKit.git
pip install -e .
```