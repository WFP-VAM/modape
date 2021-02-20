#!/bin/bash 
SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname "$SCRIPT")

if [ ! -d "$SCRIPTPATH/../storage" ]
then
    mkdir $SCRIPTPATH/../storage
fi