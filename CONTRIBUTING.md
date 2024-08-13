# EMToolKit


## Installation Instructions for Development

<!-- ### Linux

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
    ``` -->

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
By following these instructions, you should be able to set up the development environment for EMToolKit on macOS. If you have any questions or need further assistance, feel free to reach out!
