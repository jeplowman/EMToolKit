import os
import sys
import logging

# Add the project root to the Python path dynamically
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

# Basic project information
project = 'EMToolKit'
author = 'Joseph Plowman'
release = '0.1.0'

# Sphinx extensions
extensions = [
    'myst_parser',
    'sphinx.ext.autosummary',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinx_autodoc_typehints',
    'sphinx_gallery.gen_gallery',
    'nbsphinx',
    'sphinx.ext.mathjax',
    'sphinx_toggleprompt',
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

# Source suffixes for reStructuredText and Markdown
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

# Sphinx-Gallery configuration
sphinx_gallery_conf = {
    'examples_dirs': [os.path.abspath('./examples')],
    'gallery_dirs': [os.path.abspath('./.auto_examples')],
    'filename_pattern': r'.*\.(ipynb|py)$',
    'capture_repr': ('_repr_html_', '__repr__'),
    'image_scrapers': ('matplotlib',),
    'doc_module': ('EMToolKit',),
    'backreferences_dir': 'gen_modules/backreferences',
    'show_memory': True,
    'download_all_examples': True,
    'thumbnail_size': (400, 280),
    'min_reported_time': 1,
    'line_numbers': True,
    'remove_config_comments': True,
}

# Setup function to add JS files
def setup(app):
    try:
        app.add_js_file('open_external_links_in_new_tab.js')
        logging.info("Added JavaScript file for external links.")
    except Exception as e:
        logging.error(f"Error adding JS file: {e}")