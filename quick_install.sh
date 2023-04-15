#!/bin/bash

#Following command will create a virtual environment for the tool. It is highly recommended as it avoids any
#library compatibility issues. You can name anything of your choice instead of 'dmra_tool'. This is required only
#once and the name should be unique. You can only activate it for the other times when you login to your server and want to use the tool.

echo "Step 1: Creating Virtual Environment: "
conda create -n dmra_tool python==3.9.16
#conda remove -n dmra_tool --all

#Second step is to activate the created virtual environment.
echo "Step 2: Activating Virtual Environment: "
conda activate dmra_tool

#now install pip, this is also required only once
echo "Step 3: Installing Dependencies (This may take several minutes): "
#conda install pip

# next is to install package dependency, we have already provided requirments file with the package. Run following
#pip install -r requirements.txt
conda install --name dmra_tool --file dmra_tool.yml

#install the package. This command should be run in the folder where setup.py and pyproject.toml is placed.
echo "Step 4: Installing Package: "
#python setup.py install
#python -m pip install .
pip install .
