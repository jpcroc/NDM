#!/bin/env python

"""
create fr/LC_MESSAGES/iradinaGUI.mo from r/LC_MESSAGES/iradinaGUI.po
"""

import polib
txt = open('fr/LC_MESSAGES/iradinaGUI.po', 'r').read()
po = polib.pofile(pofile=txt, encoding='utf-8')
po.save_as_mofile('fr/LC_MESSAGES/iradinaGUI.mo')
print("OK translate.py seems to be done")
