#!/bin/sh
if [ -f "../../setup.py" ]
then
    # When developing within the context of WFP-VAM MODAPE fork (cwd => chain/.devcontainer):
    cd ../..
else
    # When developing/deploying with WFP-VAM MODAPE installed as dependecy:
    echo "Cloning MODAPE from WFP-VAM repository..."
    cd /var/tmp && git clone https://github.com/WFP-VAM/modape.git && cd modape && git checkout tags/v1.0 -b modape-v1.0
fi
pip install .
