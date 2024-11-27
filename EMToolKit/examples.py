import webbrowser
import subprocess
import time
import os

def example_dir():
    url = "https://emtoolkit.readthedocs.io/en/latest/examples/GALLERY_HEADER.html"
    webbrowser.open(url)


notebook_path = "EMToolKit_top_example_07252010.ipynb"


def start_jupyter_notebook(notebook_path, open_browser=True):
    # Resolve the full path of the notebook
    notebook_path = os.path.abspath(notebook_path)
    notebook_dir = os.path.dirname(notebook_path)

    # Construct the command to start Jupyter server with the specific notebook
    command = ["jupyter", "notebook", notebook_path]

    # Start the Jupyter server
    process = subprocess.Popen(command)

    # Give the server some time to start
    time.sleep(2)

    # Open the notebook URL if requested
    if open_browser:
        url = f"http://localhost:8888/notebooks/{os.path.basename(notebook_path)}"
        webbrowser.open(url)

    return process

def stop_jupyter_notebook(process):
    process.terminate()


def example_run():
    # Start the Jupyter notebook
    jupyter_process = start_jupyter_notebook(notebook_path)

    # To keep the script running and the server open
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        print("Stopping the Jupyter server...")
        stop_jupyter_notebook(jupyter_process)


if __name__ == "__main__":
    examples()
