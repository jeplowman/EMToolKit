
# Android Setup for Running EMToolKit on JupyterLab

If you're interested in improving the mobile version of the dashboard, start here!

## 1. Install F-Droid
- Navigate to [https://f-droid.org/packages/com.termux/](https://f-droid.org/packages/com.termux/).
- Install F-Droid on your Android device.

## 2. Install Termux
- Open the F-Droid app and search for Termux.
- Install Termux from F-Droid.
- Open Termux and update the package manager:
  ```sh
  pkg update && pkg upgrade
  apt update && apt upgrade
  ```

## 3. Install Pip, Git, and Other Dependencies
- Install Pip and Git using Termux:
  ```sh
  apt install python-pip git
  ```
- Install additional required packages:
  ```sh
  apt install python-dev libzmq clang libcrypt-dev freetype freetype-dev openssl
  ```
- Upgrade Pip:
  ```sh
  pip install --upgrade pip
  ```

## 4. Clone the Repository
- Clone the EMToolKit repository to your device:
  ```sh
  git clone https://github.com/jeplowman/EMToolKit.git
  ```
- Navigate into the repository directory:
  ```sh
  cd EMToolKit
  ```
- Configure Git to recognize the repository as safe:
  ```sh
  git config --global --add safe.directory $(pwd)
  ```

## 5. Install JupyterLab
- Install JupyterLab via pip:
  ```sh
  pip install jupyterlab
  ```

## 6. Set Up and Run JupyterLab
- Start JupyterLab on your mobile device:
  ```sh
  jupyter lab --ip=0.0.0.0 --no-browser
  ```
- JupyterLab will start, and you’ll see output similar to:
  ```
  http://localhost:8888/lab?token=<unique-token>
  ```

## 7. Access JupyterLab in Chrome
- Open the Chrome browser on your Android device.
- Enter the URL provided by JupyterLab in Termux, including the token (e.g., `http://localhost:8888/lab`).
- This will open JupyterLab in Chrome, allowing you to run and interact with the EMToolKit software.

## 8. Running the EMToolKit Software
- Within JupyterLab in Chrome, navigate to the directory where the EMToolKit repository was cloned (`EMToolKit`).
- Open a notebook or create a new one, and you can now start running the EMToolKit software directly on your mobile device.
- To run your software, use the notebooks or Python scripts as you would on a desktop environment.

## 9. Stopping JupyterLab
- To stop the JupyterLab server, go back to the Termux terminal and press `Ctrl + C`.

## 10. (Optional) Install the GitHub App
- For easier access to GitHub, you can also install the GitHub app from the Google Play Store to manage your repositories directly on your mobile device.
