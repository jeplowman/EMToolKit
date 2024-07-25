# EMToolKit

Repository for EMToolKit project for computing Solar Differential Emission Measures in Python.

To run the example notebook, you'll need Python set up with Sunpy, Numpy, Astropy, and the latest version of NDCube (see [NDCube Installation Guide](https://docs.sunpy.org/projects/ndcube/en/latest/installation.html)). After that, the example should work with the existing directory structure. Let me know if not, of course!

## Installation Instructions for Development

### Prerequisites

- Ensure Pip and Conda are installed.

### Linux

1. **Install Pip and Conda**:
    ```sh
    sudo apt-get install pip

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    source ~/.bashrc
    conda init
    ```

2. **Create and Activate the Conda Environment**:
    ```sh
    conda create --name EMToolKit_env
    conda activate EMToolKit_env
    ```

3. **Install the Package**:
    ```sh
    pip install -e .
    ```

### macOS

1. **Install Conda and Create a Blank Environment**:

    - Install Conda (if not already installed):
      Follow the instructions on the [Conda installation page](https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html).

2. **Clone the Repository**:
    ```sh
    git clone https://github.com/jeplowman/EMToolKit.git
    cd EMToolKit
    ```

3. **Create and Activate the Conda Environment**:
    ```sh
    conda env create -f environment.yml
    conda activate EMToolKit_env
    ```

4. **Install the Package in Editable Mode**:
    ```sh
    pip install -e .
    ```

## Usage

Once the environment is set up, you can start using the toolkit as described in the project documentation.

## Additional Information

- **Contributing**: Contributions are welcome! Please see the [CONTRIBUTING.md](CONTRIBUTING.md) file for more details.
- **License**: This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.
- **Troubleshooting**: If you encounter any issues during the setup, please refer to the project's issue tracker or documentation.

### Example Notebook

To run the example notebook, ensure your environment includes the necessary dependencies (Sunpy, Numpy, Astropy, and NDCube). Follow the installation instructions provided above. If the example does not work as expected, please let us know.

---

By following these instructions, you should be able to set up the development environment for EMToolKit on both Linux and macOS. If you have any questions or need further assistance, feel free to reach out!