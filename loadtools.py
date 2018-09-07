#!/usr/bin/env python

import os
import numpy as np

# Generic save and load

# Appends and updates values if already existing
def saveNPY(fname, **nargs):
    data = nargs
    try:
        res, fname = loadNPY(fname)
        data.update(res)
    except IOError:

        filename, file_extension = os.path.splitext(fname)
        if file_extension is '':
            file_extension = '.npy'
        fname = filename + file_extension
        print "Creating new file: {}".format(fname)
        pass
    np.save(fname, data)

# if not file_extension is given '.npy' is appended
def loadNPY(fname):
    filename, file_extension = os.path.splitext(fname)
    if file_extension is '':
        file_extension = '.npy'
    fname = filename + file_extension
    res = np.load(fname).item()
    return res, fname



# Takes and save SS Model
def saveModel(fname, A, B, C, D):
    saveNPY(fname, A=A, B=B, C=C, D=D)

# Extracts SS  Model
def getModel(fname):
    res, _ = loadNPY(fname)
    if all (k in res for k in ('A', 'B', 'C', 'D')):
        return res['A'], res['B'], res['C'], res['D']
    else:
        raise KeyError('Not all matrizes (A, B, C, D) where found in loaded dictionary')



# Takes and saves Dynamical Control System Matrizes
def saveControl(fname, Ak, Bk, Ck, Dk, Ek, u_max):
    saveNPY(fname, Ak=Ak, Bk=Bk, Ck=Ck, Dk=Dk, Ek=Ek, u_max=u_max)

# Extracts Dynamical Control System Matrizes
def getControl(fname):
    res, _ = loadNPY(fname)
    if all (k in res for k in ('Ak', 'Bk', 'Ck', 'Dk', 'Ek', 'u_max')):
        return res['Ak'], res['Bk'], res['Ck'], res['Dk'], res['Ek'], res['u_max']
    else:
        raise KeyError('Not all matrizes (Ak, Bk, Ck, Dk, Ek, u_max) where found in loaded dictionary')



# Takes and save SS Model iekf style
def saveIEKF(fname, w, s, p, b, c, d):
    iekf = {
        'w' : w,
        's' : s,
        'p' : p,
        'b' : b,
        'c' : c,
        'd' : d
    }
    saveNPY(fname, iekf_data=iekf)

# Extracts SS  Model iekf style
def getIEKF(fname):
    res, _ = loadNPY(fname)
    if 'iekf_data' in res:
        return res['iekf_data']
    else:
        raise KeyError('IEKF data not found in loaded dictionary')

def saveDelayModel(fname, A, B, C, D, delay):
    saveModel(fname, A, B, C, D)
    saveNPY(fname, delay=delay)

def getDelayModel(fname):
    A, B, C, D = getModel(fname)
    res, _ = loadNPY(fname)
    return A, B, C, D, res['delay']
