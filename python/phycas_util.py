#!/usr/bin/python
import os

def quoteDirIfNeeded(s):
    '''returns a path surrounded with double quotes to paths if it has spaces'''
    if len(s.split()) > 1:
        if len(s.split('"')) > 1:
            return '"' + '\"'.join(s.split('"')) + '"'
        return '"' + s + '"'
    return s
