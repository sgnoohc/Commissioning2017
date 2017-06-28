#!/bin/bash

usage()
{
    echo "Usage:"
    echo ""
    echo "  $(basename $0) PATHTOROOT OUTPUTPREFIX[_output.root]"
    echo ""
    echo ""
    exit
}

if [ -z $1 ]; then usage; fi
if [ -z $2 ]; then usage; fi

root -l -b -q doAll.C+\(\"$1\",\"$2\"\)
