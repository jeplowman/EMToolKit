import webbrowser
import subprocess
import time
import os
from importlib import resources
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)

# Example module for resource management
import EMToolKit


def example_dir():
    """
    Opens the example gallery in the web browser.
    """
    url = "https://emtoolkit.readthedocs.io/en/latest/examples/GALLERY_HEADER.html"
    webbrowser.open(url)
    logging.info("Opened example gallery in the web browser.")


def get_notebook_path(filename="EMToolKit_top_example_07252010.ipynb"):
    """
    Retrieves the path to a specified Jupyter Notebook located in the package.

    Args:
        filename (str): Name of the notebook file.

    Returns:
        str: Absolute path to the notebook.
    """
    try:
        notebook_path = resources.files("EMToolKit.examples").joinpath(filename)
        if not notebook_path.exists():
            raise FileNotFoundError(f"Notebook not found in package: {notebook_path}")
        return str(notebook_path)
    except FileNotFoundError as e:
        logging.error(e)
        raise


def start_jupyter_notebook(notebook_path, open_browser=True):
    """
    Launches a Jupyter Notebook using subprocess.

    Args:
        notebook_path (str): Path to the Jupyter Notebook file.
        open_browser (bool): Whether to open the notebook in the browser.

    Returns:
        subprocess.Popen: The Jupyter server process, if started successfully.
    """
    notebook_dir = os.path.dirname(notebook_path)
    if not os.path.exists(notebook_path):
        raise FileNotFoundError(f"The specified notebook does not exist: {notebook_path}")

    try:
        os.chdir(notebook_dir)
        logging.info("Launching Jupyter Notebook...")
        command = ["jupyter", "notebook", notebook_path]
        if not open_browser:
            command.append("--no-browser")
        process = subprocess.Popen(command)
        time.sleep(2)  # Allow the server time to initialize
        logging.info("Notebook launched successfully!")
        return process
    except Exception as e:
        logging.error(f"Failed to launch notebook: {e}")
        raise


def stop_jupyter_notebook(process):
    """
    Stops the Jupyter server process.

    Args:
        process (subprocess.Popen): The server process to terminate.
    """
    if process:
        logging.info("Terminating the Jupyter server...")
        process.terminate()


def example_run():
    """
    Main function to start the example notebook and keep it running.
    """
    try:
        # Get the path to the included notebook
        notebook_path = get_notebook_path()

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


if __name__ == "__main__":
    example_dir()  # Open the example directory in the browser
    example_run()  # Run the example notebook