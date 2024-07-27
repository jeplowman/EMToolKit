import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

project = 'EMToolKit'
author = 'Joseph Plowman'
release = '0.1.0'

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

import subprocess

def run_apidoc(_):
    source_dir = os.path.abspath('EMToolKit')
    output_dir = os.path.abspath('docs/source')
    subprocess.run(['sphinx-apidoc', '-o', output_dir, source_dir])

def setup(app):
    app.connect('builder-inited', run_apidoc)