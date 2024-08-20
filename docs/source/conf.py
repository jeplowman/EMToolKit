import os
import sys
from sphinx_gallery.sorting import FileNameSortKey
import logging

# Add the project root to the Python path dynamically
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

# Basic project information
project = 'EMToolKit'
author = 'Joseph Plowman'
release = '0.1.0'

# Sphinx extensions
# extensions = [
#     'myst_parser',  # If you need Markdown support
#     # 'sphinx.ext.autosummary',
#     'sphinx.ext.autodoc',
#     # 'sphinx.ext.napoleon',
#     # 'sphinx.ext.todo',
#     'sphinx.ext.viewcode',
#     # 'sphinx_autodoc_typehints',
#     # 'sphinx_gallery.gen_gallery',  # For handling Python scripts
#     'nbsphinx',  # For handling Jupyter Notebooks
#     # 'sphinx.ext.mathjax',
#     # 'sphinx_toggleprompt',
# ]
extensions = [
    'nbsphinx',  # To process Jupyter Notebooks
    'sphinx.ext.autodoc',  # To generate documentation from .py docstrings
    'sphinx.ext.viewcode',  # To add links to highlighted source code
    'sphinx_gallery.gen_gallery',  # To handle Python scripts as gallery examples
    'myst_parser',  # If you need Markdown support
]
add_module_names = False

# MyST parser extensions
myst_enable_extensions = [
    "amsmath",
    "dollarmath",
    "deflist",
]

# Autosummary settings
autosummary_generate = True
autosummary_imported_members = True

# Paths and templates
templates_path = ['_templates']
exclude_patterns = []
html_static_path = ['_static']

# HTML theme and options
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'collapse_navigation': False,
    'sticky_navigation': True,
    'titles_only': True,
}

# nbsphinx settings
nbsphinx_allow_errors = True

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
    # '.ipynb': 'nbsphinx',  # Recognize Jupyter Notebooks
    # '.py': 'python',  # Recognize Python files (note: you need an extension to process them)
}

# Sphinx-Gallery configuration for handling Python scripts
# sphinx_gallery_conf = {
#     'examples_dirs': ['./examples'],  # Path to your example scripts
#     'gallery_dirs': ['../build/auto_examples'],  # Where to save the gallery output
#     'filename_pattern': r'.*\.(ipynb|py)$',  # Match both notebooks and .py files
#     'capture_repr': ('_repr_html_', '__repr__'),
#     'image_scrapers': ('matplotlib',),
#     'doc_module': ('EMToolKit',),
#     'backreferences_dir': '../build/gen_modules/backreferences',
#     'show_memory': True,
#     'download_all_examples': True,
#     'thumbnail_size': (400, 280),
#     'min_reported_time': 1,
#     'line_numbers': True,
#     'remove_config_comments': True,
#     'within_subsection_order': FileNameSortKey,
# }

sphinx_gallery_conf = {
    'examples_dirs': './examples',  # Path to example scripts
    'gallery_dirs': '../build/auto_examples',  # Output directory for the gallery
    'filename_pattern': r'.*\.py$',  # Only match Python scripts
}


# Setup function to add JS files
def setup(app):
    try:
        app.add_js_file('open_external_links_in_new_tab.js')
        logging.info("Added JavaScript file for external links.")
    except Exception as e:
        logging.error(f"Error adding JS file: {e}")