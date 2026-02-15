#!/bin/bash

if [ -f "super_secret_key.txt" ]; then
    echo "Found it!" > found.log
else
    echo "I did my best"
    ls super_secret_key.txt 2>/dev/null
fi
