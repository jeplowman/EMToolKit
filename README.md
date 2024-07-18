# EMToolKit
Repository for EMToolKit project for computing Solar Differential Emission Measures in Python

To run the example notebook, you'll need python set up with Sunpy, Numpy, Astropy, and the latest version of NDCube (see https://docs.sunpy.org/projects/ndcube/en/latest/installation.html). After that, example should work with the existing directory structure. Let me know if not, of course!


## Installation Instructions for Development
### Linux
Ensure Pip and Conda are installed

    $ apt-get install pip

    $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

    $ bash Miniconda3-latest-Linux-x86_64.sh

    $ source ~./bashrc

    $ conda init

Create the environment

    conda create --name EMToolKit_env
    conda activate EMToolKit_env

Install the Package

    pip install

### Macos

# Install conda and make a blank environment
conda init #If you've never done it before, initialize the shell
conda create --name emtoolkitenvarm # Make your new environment
conda activate emtoolkitenvarm # Activate the environment
conda install python
git clone https://github.com/jeplowman/EMToolKit.git


pip install -e .


