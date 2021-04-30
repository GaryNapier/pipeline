#! /bin/bash

if [ ! -d ~/bin/ ]; then
    mkdir ~/bin/
fi

cp ./shell_scripts/* ~/bin/
cp ./r_scripts/* ~/bin/
cp ./python_scripts/* ~/bin/

chmod +x ~/bin/*.R
chmod +x ~/bin/*.sh
chmod +x ~/bin/*.py
