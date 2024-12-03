import os
import sys
import argparse
import json
import time
import logging
import webbrowser
import subprocess
from pathlib import Path
import importlib.resources as pkg_resources
from .examples.examples import example_dir, example_run

def start_jupyter_notebook(notebook_path):
    """
    Start the Jupyter notebook.

    Args:
        notebook_path (str): Path to the Jupyter notebook to open.

    Returns:
        Popen: The running Jupyter process.
    """
    return subprocess.Popen(["jupyter", "notebook", notebook_path])

def stop_jupyter_notebook(jupyter_process):
    """
    Stop the Jupyter notebook process.

    Args:
        jupyter_process (Popen): The running Jupyter process.
    """
    jupyter_process.terminate()

def get_notebook_path(target_directory, notebook_name):
    """
    Get the path to the newly created notebook.

    Args:
        target_directory (str): Directory where the notebook is created.
        notebook_name (str): Name of the notebook file.

    Returns:
        str: The path to the notebook.
    """
    return os.path.join(target_directory, notebook_name)

def open_in_file_browser(directory):
    """
    Open the specified directory in the system's file browser.

    Args:
        directory (str): The directory to open.
    """
    if os.name == 'nt':  # Windows
        os.startfile(directory)
    elif os.name == 'posix':  # macOS and Linux
        subprocess.Popen(['open', directory])  # 'open' is macOS specific
        # You might need 'xdg-open' for Linux: subprocess.Popen(['xdg-open', directory])

def example_run(notebook_path, no_preview=False):
    """
    Main function to start the example notebook, open in browser, and keep it running.

    Args:
        notebook_path (str): Path to the notebook.
        no_preview (bool): If True, do not open the notebook in the web browser.
    """
    try:
        notebook_dir = str(Path(notebook_path).parent)

        # Open the directory in the file browser
        open_in_file_browser(notebook_dir)

        if not no_preview:
            # Open the notebook in the default web browser
            webbrowser.open_new(f'file://{os.path.abspath(notebook_path)}')

        # Start the Jupyter notebook
        jupyter_process = start_jupyter_notebook(notebook_path)

        # Keep the script running and the server open
        logging.info("Press Ctrl+C to stop the Jupyter server.")
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        logging.info("Stopping the Jupyter server...")
        if 'jupyter_process' in locals():
            stop_jupyter_notebook(jupyter_process)
    except Exception as e:
        logging.error(f"Error encountered: {e}")

def create_new_notebook(target_directory, *, images_directory=None, notebook_name=None):
    """
    Creates a new Jupyter Notebook in the specified target directory.

    Args:
        target_directory (str): Directory to save the new notebook.
        images_directory (str, optional): Path to the directory containing images. Defaults to None.
        notebook_name (str, optional): Name of the notebook file. Defaults to "new_notebook.ipynb".
    """
    try:
        with pkg_resources.open_text(__package__, 'notebook_template.json') as file:
            notebook_content = json.load(file)
    except FileNotFoundError:
        print(f"Template file not found.")
        return None
    except json.JSONDecodeError:
        print(f"Error decoding JSON from template file.")
        return None
    except Exception as e:
        print(f"Failed to load template: {e}")
        return None

    # Ensure the target directory exists
    os.makedirs(target_directory, exist_ok=True)

    # Construct the notebook file path
    if notebook_name is None:
        notebook_name = f"{os.path.basename(os.getcwd())}_NB.ipynb"
    notebook_path = os.path.join(target_directory, notebook_name)

    # Ensure the path does not conflict with an existing directory
    if os.path.isdir(notebook_path):
        raise IsADirectoryError(f"Cannot create notebook: {notebook_path} is a directory.")

    try:
        with open(notebook_path, "w") as notebook_file:
            json.dump(notebook_content, notebook_file, indent=4)
        print(f"Notebook successfully created at {notebook_path}")
        return notebook_path
    except Exception as e:
        print(f"Failed to create notebook: {e}")
        return None

def main():
    """
    Entry point for the EMToolKit CLI tool.
    """
    parser = argparse.ArgumentParser(
        description="EMToolKit CLI - Manage notebooks and explore directories."
    )
    subparsers = parser.add_subparsers(
        dest="command",
        help="Available commands",
        title="Primary Commands",
    )

    # Primary command: create
    create_parser = subparsers.add_parser(
        "create",
        help="Create a new notebook in the target directory (default command)",
    )
    create_parser.add_argument(
        "target_directory",
        nargs="?",
        type=str,
        help="Directory where the notebook will be created (defaults to the current directory)",
    )
    create_parser.add_argument(
        "images_directory",
        nargs="?",
        type=str,
        help="Path to the images directory (defaults to the target directory)",
    )
    create_parser.add_argument(
        "--no-preview",
        action="store_true",
        help="Do not open the notebook in a web browser",
    )

    # Additional commands grouped under a subsection
    parser.add_argument(
        "--examples",
        action="store_true",
        help="Show the example directory in the readthedocs",
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Run the example notebook in this terminal",
    )

    # Parse the arguments
    args = parser.parse_args()

    # Handle the create command
    if args.command == "create" or not args.command:
        target_directory = args.target_directory or os.getcwd()
        images_directory = args.images_directory or target_directory
        notebook_path = create_new_notebook(target_directory, images_directory=images_directory)
        if not notebook_path:
            sys.exit(1)
        else:
            example_run(notebook_path, no_preview=args.no_preview)
    # Handle additional commands
    elif args.examples:
        example_dir()
    elif args.test:
        example_run(None, no_preview=True)  # Assuming example_run handles its own notebook path
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
