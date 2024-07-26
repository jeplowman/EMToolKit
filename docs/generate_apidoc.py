# docs/generate_apidoc.py

import os
import subprocess

# Define the source and output directories
source_dir = 'EMToolKit'
output_dir = 'docs/source'

# Run sphinx-apidoc
subprocess.run(['sphinx-apidoc', '-o', output_dir, source_dir])