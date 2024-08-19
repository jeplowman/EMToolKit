
# EMToolKit Development Guide

By following these instructions, you should be able to set up the development environment for EMToolKit on your preferred operating system. If you have any questions or need further assistance, feel free to reach out.

## macOS
This is the most stable development platform currently.

1. **Install Conda**:
   - If Conda is not already installed, follow the instructions on the [Conda installation page](https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html).

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

---

## Android
If you're interested in improving the mobile version of the dashboard, start here!

1. **Install F-Droid**:
   - Navigate to [https://f-droid.org/packages/com.termux/](https://f-droid.org/packages/com.termux/).
   - Install F-Droid.

2. **Install Termux**:
   - Install Termux from F-Droid.
   - Update the package manager:
     ```sh
     pkg update && pkg upgrade
     apt update && apt upgrade
     ```

3. **Install Pip and Git**:
   - Install Pip and Git using Termux:
     ```sh
     apt install pip git
     ```
   - (Optional) Install the GitHub app from the app store.

4. **Clone the Repository**:
   - Clone the EMToolKit repository:
     ```sh
     git clone https://github.com/jeplowman/EMToolKit.git
     ```
   - Configure Git to recognize the repository as safe:
     ```sh
     git config --global --add safe.directory [path]
     ```

---
<!--
## Linux
### Untested ->

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

3. **Install the Package in Editable Mode**:
   ```sh
   pip install -e .
   ``` -->
