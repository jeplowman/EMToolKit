
# EMToolKit Development Guide

## Dev Installation
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
Check out the [Android Setup](Android_EMToolKit_Setup.md) page for more information about developing for the mobile version of the app.
