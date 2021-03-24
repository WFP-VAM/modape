#!/bin/sh
if [ -f "../../setup.py" ]
then
    cd ../..
else
    echo "Cloning MODAPE from WFP-VAM repository..."
    cd /workspaces && git clone https://github.com/WFP-VAM/modape.git && cd modape && git checkout -b v1.0rc origin/v1.0rc
fi
pip install .
