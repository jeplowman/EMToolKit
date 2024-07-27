import os
import sys
import subprocess

# Add the project root to the Python path
sys.path.insert(0, os.path.abspath('../../'))

project = 'EMToolKit'
author = 'Joseph Plowman'
release = '0.1.0'

# List of Sphinx extensions
extensions = [
    'sphinx.ext.autosummary',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx_autodoc_typehints',
]

autosummary_generate = True  # Turn on autosummary
autosummary_imported_members = True

templates_path = ['_templates']
exclude_patterns = []

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

def run_apidoc(_):
    """Generate .rst files for the Sphinx documentation."""
    try:
        source_dir = os.path.abspath('../../EMToolKit')  # Adjust the path as necessary
        output_dir = os.path.abspath('docs/source')
        subprocess.check_call(['sphinx-apidoc', '-o', output_dir, source_dir])
    except subprocess.CalledProcessError as e:
        print("sphinx-apidoc failed with exit code", e.returncode)

def setup(app):
    app.connect('builder-inited', run_apidoc)