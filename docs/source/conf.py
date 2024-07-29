import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

project = 'EMToolKit'
author = 'Joseph Plowman'
release = '0.1.0'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx_autodoc_typehints',
    'nbsphinx',
    'sphinx.ext.mathjax',
    'sphinx_gallery.gen_gallery',
]

templates_path = ['_templates']
exclude_patterns = []

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Configure sphinx-gallery
sphinx_gallery_conf = {
    'examples_dirs': 'examples',  # path to your example scripts
    'gallery_dirs': 'examples/.auto_examples',  # path to where to save generated output
}

def run_apidoc(_):
    """Generate .rst files for the Sphinx documentation."""
    try:
        source_dir = os.path.abspath('../../EMToolKit/EMToolKit')  # Adjust the path as necessary
        output_dir = os.path.abspath('docs/source')
        subprocess.check_call(['sphinx-apidoc', '-o', output_dir, source_dir])
    except subprocess.CalledProcessError as e:
        print("sphinx-apidoc failed with exit code", e.returncode)
        try:
            source_dir = os.path.abspath('../EMToolKit/EMToolKit')  # Adjust the path as necessary
            output_dir = os.path.abspath('docs/source')
            subprocess.check_call(['sphinx-apidoc', '-o', output_dir, source_dir])
        except Exception as e:
            raise e

def setup(app):
    app.connect('builder-inited', run_apidoc)