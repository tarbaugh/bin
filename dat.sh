#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "No arguments supplied"
    return	
fi


python ~/bin/strip.py "$@"
