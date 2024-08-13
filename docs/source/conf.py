import os
import sys
import subprocess

# Add the project root to the Python path
sys.path.insert(0, os.path.abspath('../../'))
print(sys.path)
print()

project = 'EMToolKit'
author = 'Joseph Plowman'
release = '0.1.0'

add_module_names = False

# List of Sphinx extensions
extensions = [
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

autosummary_generate = True  # Turn on autosummary
autosummary_imported_members = True

templates_path = ['_templates']
exclude_patterns = []

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

html_theme_options = {
    'collapse_navigation': False,
    'sticky_navigation': True,
    'titles_only': True,  # This makes the sidebar only show titles
}

nbsphinx_allow_errors = True

# Configure sphinx-gallery
sphinx_gallery_conf = {
    'examples_dirs': [os.path.abspath('./examples')],
    'gallery_dirs': [os.path.abspath('./.auto_examples')],
    'filename_pattern': r'.*\.(ipynb|py)$',  # Match both .py and .ipynb files
    'capture_repr': ('_repr_html_', '__repr__'),  # Capture HTML representations
    'image_scrapers': ('matplotlib',),  # Specify scrapers for image output
    'doc_module': ('EMToolKit',),  # Document the EMToolKit module
    'backreferences_dir': 'gen_modules/backreferences',  # Path for backreferences
    'show_memory': True,  # Show memory usage
    'download_all_examples': True,  # Allow downloading all examples in a zip file
    'thumbnail_size': (400, 280),  # Set the size of thumbnails
    'min_reported_time': 1,  # Report all examples, regardless of runtime
    'line_numbers': True,  # Show line numbers in code blocks
    'remove_config_comments': True,  # Remove comments from config files
}

def run_apidoc(app):
    """Generate .rst files for the Sphinx documentation and build HTML output."""
    if os.environ.get('SPHINX_APIDOC_RUNNING') == '1':
        print("Skipping sphinx-apidoc because it is already running.")
        return

    if 'skip-apidoc' in app.config.sphinx_run_options:
        print("Skipping sphinx-apidoc because skip-apidoc flag was set.")
        return

    try:
        os.environ['SPHINX_APIDOC_RUNNING'] = '1'

        source_dir = os.path.abspath('EMToolKit/')
        output_dir = os.path.abspath('docs/source/')
        html_dir = os.path.abspath("docs/build/html")

        subprocess.check_call(['sphinx-apidoc', '-o', output_dir, source_dir])
        subprocess.check_call(['sphinx-build', '-b', 'html', output_dir, html_dir])

    except subprocess.CalledProcessError as e:
        print(f"sphinx-apidoc or sphinx-build failed with exit code {e.returncode}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
    finally:
        os.environ.pop('SPHINX_APIDOC_RUNNING', None)

def setup(app):
    app.add_config_value('sphinx_run_options', [], 'env')
    app.connect('builder-inited', run_apidoc)