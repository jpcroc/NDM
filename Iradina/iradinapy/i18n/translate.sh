#!/bin/bash

# this script creates the binary resources file from the strings file

I18HOME=`dirname $0`

# french
echo "build French ressources, create iradinaGUI.mo"
msgfmt ${I18HOME}/fr/LC_MESSAGES/iradinaGUI.po -o ${I18HOME}/fr/LC_MESSAGES/iradinaGUI.mo
 
