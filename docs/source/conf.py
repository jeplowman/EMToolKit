import os
import sys
sys.path.insert(0, os.path.abspath('../../'))
print(sys.path)
print()

project = 'EMToolKit'
author = 'Joseph Plowman'
release = '0.1.0'

add_module_names = False


# List of Sphinx extensions
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx_autodoc_typehints',
    'nbsphinx',
    'sphinx.ext.mathjax',
    'sphinx_gallery.gen_gallery',
    'sphinx_gallery.load_style'
]

templates_path = ['_templates']
exclude_patterns = []

html_theme = 'sphinx_rtd_theme'
# html_theme = 'sphinxdoc'
# html_theme = 'furo'
html_static_path = ['_static']

html_theme_options = {
    'collapse_navigation': False,
    'sticky_navigation': True,
    'titles_only': True,  # This makes the sidebar only show titles
}

nbsphinx_allow_errors = True

# Configure sphinx-gallery
sphinx_gallery_conf = {
    'examples_dirs': ['examples'],  # paths to your notebooks and example scripts
    # 'gallery_dirs': ['../build/auto_examples'],  # paths to where to save generated output
    'gallery_dirs': ['../build/html/examples'],  # paths to where to save generated output
    # 'gallery_dirs': ['examples'],  # paths to where to save generated output
    'filename_pattern': r'.*\.(ipynb|py)$',  # Match both .py and .ipynb files
    'within_subsection_order': 'FileNameSortKey',  # Order examples by file name
    # 'default_thumb_file': 'path/to/default_thumbnail.png',  # Optional: set a default thumbnail
    'capture_repr': ('_repr_html_', '__repr__'),  # Capture HTML representations
    'image_scrapers': ('matplotlib',),  # Specify scrapers for image output
    'doc_module': ('EMToolKit',),  # Document the EMToolKit module
    'reference_url': {
        'EMToolKit': "https://github.com/jeplowman/EMToolKit/tree/main",  # Set reference URL for EMToolKit
    },
    'backreferences_dir': 'gen_modules/backreferences',  # Path for backreferences
    'show_memory': True,  # Show memory usage
    # 'junit': True,  # Enable JUnit XML output
    # 'junit_dir': 'gen_modules/junit',  # Directory for JUnit XML files
    'download_all_examples': True,  # Allow downloading all examples in a zip file
    'thumbnail_size': (400, 280),  # Set the size of thumbnails
    'min_reported_time': 0,  # Report all examples, regardless of runtime
    # 'binder': {
    #     'org': 'your-github-org',  # GitHub organization or user
    #     'repo': 'your-repo',  # GitHub repository
    #     'branch': 'main',  # Branch
    #     'binderhub_url': 'https://mybinder.org',  # BinderHub URL
    #     'dependencies': ['./requirements.txt'],  # Path to requirements file
    #     'notebooks_dir': 'notebooks',  # Path to notebooks directory
    #     'use_jupyter_lab': True,  # Use JupyterLab interface
    # },
    'compress_images': ('images', 'thumbnails'),  # Compress images
    'line_numbers': True,  # Show line numbers in code blocks
    'remove_config_comments': True,  # Remove comments from config files
}

import os, subprocess
print(f"Sphinx root directory: {os.path.abspath(os.path.dirname(__file__))}\n")

def run_apidoc(_):
    """Generate .rst files for the Sphinx documentation."""
    try:
        source_dir = os.path.abspath('EMToolKit/')  # Adjust the path as necessary
        output_dir = os.path.abspath('docs/source/')
        subprocess.check_call(['sphinx-apidoc', '-o', output_dir, source_dir])
    except subprocess.CalledProcessError as e:
        print("sphinx-apidoc failed with exit code", e.returncode)
        try:
            source_dir = os.path.abspath('EMToolKit/EMToolKit')  # Adjust the path as necessary
            output_dir = os.path.abspath('docs/source')
            subprocess.check_call(['sphinx-apidoc', '-o', output_dir, source_dir])
        except Exception as e:
            raise e

def setup(app):
    app.connect('builder-inited', run_apidoc)