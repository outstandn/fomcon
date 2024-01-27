# FOMCONpy
FOMCONpy is a new fractional-order modelling and control toolbox for Python. It is an extension of the existing FOMCON toolbox for MATLAB, but this time aiming at Python users and the Internet of Things (IoT) community. Just like the original toolbox, it offers a set of tools for researchers in the field of fractional-order control.

This toolbox is dependent on control 0.8.1 library for python (https://python-control.readthedocs.io/en/0.8.1/intro.html). Other dependencies include pyqt5, numpy, scipy, pandas, xlrd and slycot (https://github.com/python-control/Slycot).
if you are using a windows OS you might have challenges installing slycot because there will be need to install a fotran compiler to build slycot from the repository.

recommended fotran compiler is Intel Math Kernel Library 2019.0.117, it works perfectly with anaconda 2019.3

Python environment used for this development is anaconda 2019.3, but pycharm should word perfectly. The "environment.uml" file contains the python library packages used in the development of this software. 

# Setup guide for Windows
- item Install Anaconda. (https://repo.anaconda.com/archive/Anaconda3-2019.03-Windows-x86_64.exe was used during development and ensure that python is added to ssytems path during installation)
- verify conda and python is install properly by running the commands below in a terminal
  - conda --version
  - python --version
- Install Git.
- Navigate using the command prompt to any desired directory of choice
- Clone FOMCONpy from its repository (git clone https://github.com/outstandn/fomcon.git).
- Change current working directory to "fomconpy": cd fomcon
- Create "fomcon" environment:
  - conda create -n fomcon python=3.7.7 numpy=1.18.1 scipy=1.4.1 pandas=1.0.3 xlrd=1.2.0 matplotlib=3.1.3 pip=20.0.2 pyqt=5.9.2
- Activate 'fomcon' environemnt
  - conda activate fomcon
- install other needed libraries in same 'fomcon' environment using command below
  -pip install control addict

# Setup guide for Linux
Open a terminal and enter the following commands:
- sudo apt-get install git python3-pyqt5 python3-numpy python3-scipy python3-pandas python3-xlrd python3-matplotlib python3-pip python3-select
- pip3 install control addict
- Navigate using the command prompt to any desired directory of choice
- Clone FOMCONpy from its repository(git clone https://github.com/outstandn/fomcon.git).
- Change current working directory to "fomconpy": cd fomcon

# Using Fomconpy's GUI
The FOMCONpy has three modules

## Fractional-Order System Analysis Module
The fractional-order system viewer can be accessed by starting a command prompt, navigating to the directory of your downloaded source, then run the follwoing command:
- python pyfomcon.py

## Fractional-Order System Identificaiton Module
The identification tool can be accessed by starting a command prompt, navigate to the directory of your downloaded FOMCONpy source code, then run the following commands:
- python pyfotfid.py

## Fractional-Order System Control Module
The tool can be accessed by starting two command prompt, navigate to the directory of your downloaded FOMCONpy source code, then run the following commands:
- python pyfopidopt.py

Open another command prompt, navigate to the directory of your downloaded FOMCONpy source code, then run the following command to start teh Destributed Control Server: 
- python controlServer.py
