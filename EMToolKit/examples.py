import webbrowser
import subprocess
import time
import os

notebook_path = "EMToolKit_top_example_07252010.ipynb"

def example_dir():
    """
    Opens the example gallery in the web browser.
    """
    url = "https://emtoolkit.readthedocs.io/en/latest/examples/GALLERY_HEADER.html"
    webbrowser.open(url)

def start_jupyter_notebook(notebook_path, open_browser=True):
    """
    Launches a Jupyter Notebook using subprocess.

    Args:
        notebook_path (str): Path to the Jupyter Notebook file.
        open_browser (bool): Whether to open the notebook in the browser.

    Returns:
        subprocess.Popen: The Jupyter server process, if started successfully.
    """
    notebook_path = os.path.abspath(notebook_path)
    notebook_dir = os.path.dirname(notebook_path)

    if not os.path.exists(notebook_path):
        raise FileNotFoundError(f"The specified notebook does not exist: {notebook_path}")

    try:
        os.chdir(notebook_dir)
        print("Attempting to launch notebook using subprocess...")
        command = ["jupyter", "notebook", notebook_path]
        if not open_browser:
            command.append("--no-browser")
        process = subprocess.Popen(command)
        time.sleep(2)  # Allow the server time to initialize
        print("Notebook launched successfully!")
        return process
    except Exception as e:
        print(f"Failed to launch notebook: {e}")
        raise

def stop_jupyter_notebook(process):
    """
    Stops the Jupyter server process.

    Args:
        process (subprocess.Popen): The server process to terminate.
    """
    if process:
        print("Terminating the Jupyter server...")
        process.terminate()

def example_run():
    """
    Main function to start the example notebook and keep it running.
    """
    try:
        # Start the Jupyter notebook
        jupyter_process = start_jupyter_notebook(notebook_path)

        # Keep the script running and the server open
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        print("Stopping the Jupyter server...")
        if jupyter_process:
            stop_jupyter_notebook(jupyter_process)
    except Exception as e:
        print(f"Error encountered: {e}")

if __name__ == "__main__":
    example_dir()  # Open the example directory in the browser
    example_run()  # Run the example notebook