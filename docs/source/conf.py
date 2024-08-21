import os
import sys
import logging
from sphinx_gallery.sorting import FileNameSortKey

# Add the project root to the Python path dynamically
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

# -------------------------------------------------------------------------
# Basic Project Information
# -------------------------------------------------------------------------
project = 'EMToolKit'
author = 'Joseph Plowman'
release = '0.1.0'

# -------------------------------------------------------------------------
# Sphinx Extensions
# -------------------------------------------------------------------------
extensions = [
    'nbsphinx',  # To process Jupyter Notebooks
    'sphinx.ext.autodoc',  # To generate documentation from .py docstrings
    'sphinx.ext.viewcode',  # To add links to highlighted source code
    'sphinx_gallery.gen_gallery',  # To handle Python scripts as gallery examples
    'myst_parser',  # Markdown support
    'sphinx.ext.autosummary',  # Automatic API summaries
    'sphinx.ext.napoleon',  # Support for Google/NumPy style docstrings
    'sphinx.ext.todo',  # TODO directives
    'sphinx_autodoc_typehints',  # Include type hints in documentation
    'sphinx.ext.mathjax',  # Render LaTeX equations
    'sphinx_toggleprompt',  # Toggle prompts in code blocks
    'sphinx_copybutton',  # Copy button in code blocks
    'sphinx.ext.githubpages',  # Ensure GitHub Pages works correctly
]

# -------------------------------------------------------------------------
# MyST Parser Settings (for Markdown support)
# -------------------------------------------------------------------------
myst_enable_extensions = [
    "amsmath",  # LaTeX math in Markdown
    "dollarmath",  # Use $ for inline math
    "deflist",  # Definition lists in Markdown
]

# -------------------------------------------------------------------------
# Autosummary Configuration
# -------------------------------------------------------------------------
autosummary_generate = True  # Automatically generate API summaries
autosummary_imported_members = True  # Include imported members in the summaries

# -------------------------------------------------------------------------
# TODO Configuration
# -------------------------------------------------------------------------
todo_include_todos = True  # Include TODOs in the generated documentation

# -------------------------------------------------------------------------
# Paths and Templates
# -------------------------------------------------------------------------
templates_path = ['_templates']  # Path to custom templates
exclude_patterns = []  # Patterns to exclude from the build
html_static_path = ['_static']  # Path to static files

# -------------------------------------------------------------------------
# HTML Theme and Options
# -------------------------------------------------------------------------
html_theme = 'sphinx_rtd_theme'  # Use the Read the Docs theme
html_theme_options = {
    'collapse_navigation': False,  # Do not collapse navigation
    'sticky_navigation': True,  # Keep the navigation sticky
    'titles_only': True,  # Only show titles in the navigation bar
}

# -------------------------------------------------------------------------
# nbsphinx Configuration (for Jupyter Notebooks)
# -------------------------------------------------------------------------
nbsphinx_allow_errors = True  # Allow errors in notebook execution
nbsphinx_execute = 'never'  # execute notebooks during the build

# -------------------------------------------------------------------------
# Source Suffixes
# -------------------------------------------------------------------------
source_suffix = {
    '.rst': 'restructuredtext',  # reStructuredText files
    '.md': 'markdown',  # Markdown files
    # '.ipynb': 'nbsphinx',  # Uncomment to explicitly recognize Jupyter Notebooks
    # '.py': 'python',  # Uncomment if processing Python files directly (not recommended)
}

# -------------------------------------------------------------------------
# Sphinx-Gallery Configuration (for Python scripts as examples)
# -------------------------------------------------------------------------
sphinx_gallery_conf = {
    'examples_dirs': './examples',  # Path to example scripts
    'gallery_dirs': '../build/auto_examples',  # Output directory for the gallery
    'filename_pattern': r'.*\.py$',  # Only match Python scripts
    # Uncomment the following if you need to handle .ipynb files as well
    # 'filename_pattern': r'.*\.(ipynb|py)$',
    # Additional configuration options (e.g., thumbnail size, line numbers)
    # 'capture_repr': ('_repr_html_', '__repr__'),
    # 'image_scrapers': ('matplotlib',),
    # 'doc_module': ('EMToolKit',),
    # 'backreferences_dir': '../build/gen_modules/backreferences',
    # 'show_memory': True,
    # 'download_all_examples': True,
    # 'thumbnail_size': (400, 280),
    # 'min_reported_time': 1,
    # 'line_numbers': True,
    # 'remove_config_comments': True,
    # 'within_subsection_order': FileNameSortKey,
}

# -------------------------------------------------------------------------
# Setup Function to Add Custom JavaScript Files
# -------------------------------------------------------------------------
def setup(app):
    try:
        app.add_js_file('open_external_links_in_new_tab.js')
        logging.info("Added JavaScript file for external links.")
    except Exception as e:
        logging.error(f"Error adding JS file: {e}")