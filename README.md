# VMI-RAE MCP coincidense Data Analysis

Description

## Installation

The program is currently only available for Linux.

### Dependacies

In order to install the program following dependencies must be satisfied:

[Cmake](https://cmake.org/download/)

[ROOT](https://root.cern/install/)

[Python3](https://www.python.org/downloads/)

Latex (any will do)

ROOT must be compiled with C++17.

One liner for Cmake, Python3, Latex and ROOT dependancies
```bash
sudo apt-get -y install python3 build-essential libssl-dev cmake texmaker dpkg-dev g++ gcc binutils libx11-dev libxpm-dev libxft-dev libxext-dev
```

ROOT can be installed in various ways. The following is an example, building it from source, setting C++17 as standard.

```bash
git clone --branch latest-stable --depth=1 https://github.com/root-project/root.git root_src
sudo mkdir root_build && /opt/root_install cd root_build
cmake -DCMAKE_INSTALL_PREFIX=/opt/root -DCMAKE_CXX_STANDARD=17 ../root_src
sudo cmake --build . --install -j4
sudo source /opt/root/bin/thisroot.sh
```
thisroot.sh must be sourced for all instances of bash, so it is advantageous to add it to the ~.bashrc or similar.

### Setup and Run

Run the setup.py to compile and setup the program
```bash
sudo python3 setup.py
```

Now everything is done! Enjoy and run the main.py.
```bash
python3 main.py
```

## Usage

