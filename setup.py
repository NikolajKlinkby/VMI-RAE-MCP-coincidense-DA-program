#!/usr/bin/env python

import os
import sys
import subprocess
import pkg_resources
from pkg_resources import DistributionNotFound, VersionConflict
from packaging import version


def get_cmake_version():
    output = subprocess.check_output(['cmake', '--version']).decode('utf-8')
    line = output.splitlines()[0]
    version = line.split()[2]
    return (version)

def should_install_requirement(requirement):
    should_install = False
    try:
        pkg_resources.require(requirement)
    except (DistributionNotFound, VersionConflict):
        should_install = True
    return should_install


def install_packages(requirement_list):
    try:
        requirements = [
            requirement
            for requirement in requirement_list
            if should_install_requirement(requirement)
        ]
        if len(requirements) > 0:
            subprocess.check_call([sys.executable, "-m", "pip", "install", *requirements])
        else:
            print("Requirements already satisfied.")

    except Exception as e:
        print(e)

current_dir = os.path.dirname(os.path.abspath(__file__))


if (sys.version_info < (3, 0, 0, 'final', 0)):
    raise SystemExit('Python 3.x is required!')

# Required dependencies for main.py 	

if os.name == 'nt':
    dependencies = ['tk', 'Pillow', 'numpy', 'thread6', 'matplotlib', 'imageio', 'scipy']
    install_packages(dependencies)
else:
    subprocess.run(['apt-get', 'install', 'python3-tk', 'dvipng', 'cm-super'])
    dependencies = ['Pillow>=8.3.2', 'numpy', 'thread6', 'matplotlib', 'imageio', 'scipy']
    install_packages(dependencies)

# CMake required
#if version.parse(get_cmake_version()) < version.parse("3.15.0"):
#    raise SystemExit('Cmake version 3.15.0 or newer is required!')

# Required path for root.cmake
import tkinter as tk
from tkinter import filedialog
root_dir = current_dir


def confirm_close(root_, label_):
    global root_dir
    root_dir = label_['text']
    root_.quit()


def search_for_dir_path(root_, title):
    global current_dir
    file = filedialog.askdirectory(parent=root_, initialdir=current_dir, title=title)
    return file


def search_for_root_location(root_, label_):
    directory = search_for_dir_path(root_, 'Please select /path/to/root/installation')
    label_.config(text=directory)


root = tk.Tk()
label = tk.Label(root, text='Root location: ')
label.grid(row=0, column=0, sticky=tk.E)
dir_label = tk.Label(root, text=current_dir)
dir_label.grid(row=0, column=1, sticky=tk.W)
dir_button = tk.Button(root, text="Choose directory", command=lambda: search_for_root_location(root, dir_label))
dir_button.grid(row=0, column=2, sticky=tk.E)

button = tk.Button(root, text="Click to continue",
                   command=lambda: confirm_close(root, dir_label))
button.grid(row=1, column=1)
root.mainloop()

if not os.path.exists(root_dir+"/ROOTConfig.cmake"):
    raise SystemError('No ROOTConfig.cmake in chosen folder')

# TODO make this windows friendly
def compile_cpp_file(folder):
    global current_dir
    global root_dir
    build = "/Source"+folder+"/build"
    if not os.path.exists(current_dir + build):
        os.mkdir(current_dir + build)
    call = 'cmake -B' + current_dir.replace(" ", "\ ") + build.replace(" ", "\ ") + ' -S' + current_dir.replace(" ","\ ") + build[:-6].replace(" ", "\ ") + ' -DROOT_DIR=' + root_dir.replace(" ", "\ ")
    os.system(call)
    call = 'cmake --build ' + current_dir.replace(" ", "\ ") + build.replace(" ", "\ ")
    os.system(call)


try:
    sources = ["/Load data", "/Event Trig data", "/VMI Boot", "/VMI filter", "/VMI GIF"]

    for src in sources:
        compile_cpp_file(src)

except:
    raise SystemError('Could not compile cpp files with cmake')
