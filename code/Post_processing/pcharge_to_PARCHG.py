#!/usr/bin/env python
# coding: utf-8

import os
import sys
current_path = os.path.dirname(os.path.abspath(__file__))
parent_path = os.path.dirname(current_path)
sys.path.append(parent_path)
import functions.aseFunction as af

data = af.pchargeLoad()

for i in data.chgInfo.keys():
    filename = "PARCHG.{}.ALLK".format("%04d" % data.bandIndex[i])
    af.write_CHGCAR(structure = data.atoms, grid = data.grid, chgInfo = data.chgInfo[i], filename = filename)