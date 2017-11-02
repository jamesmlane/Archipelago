# ----------------------------------------------------------------------------
#
# TITLE - misc.py
# AUTHOR - James Lane
# PROJECT - archipelago
# CONTENTS:
#	1.
#
# ----------------------------------------------------------------------------
#
# Docstrings and metadata:
'''
Miscellaneous functions for the Archipelago project.
'''

__author__ = "James Lane"

#Imports
import re

def sort_nicely(l):
    '''
    sort_nicely:

    Sort a bunch of files as a human would rather than how a dumb computer
    would

    Args:
        l (str array) - array of filenames to sort

    Returns:
        s - sorted array of
    '''
    def tryint(s):
            try:
                    return int(s)
            except:
                    return s
    #def
    def alphanum_key(s):
            return [ tryint(c) for c in re.split('([0-9]+)', s) ]
    #def
    l.sort(key=alphanum_key)
