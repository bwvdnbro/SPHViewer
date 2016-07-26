#! /usr/bin/python

################################################################################
# This file is part of SPHViewer
# Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#
# SPHViewer is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SPHViewer is distributed in the hope that it will be useful,
# but WITOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with SPHViewer. If not, see <http://www.gnu.org/licenses/>.
################################################################################

import numpy as np
import scipy as sp
import chyplot
from enums import T_gas, T_star
import matplotlib.pyplot as plt
import hyplot.comp.PGlobals as glob
import struct
import sys

reader = chyplot.CDataGadget()
print "reading snapshot", sys.argv[2][-4:]
reader.setFilename(sys.argv[2])
# read the file into a CDataBlock
data = reader.readFile()
#data.rcom(True, T_star, 0, 0, 10)
#data.vcom(True, T_star)

ngas = data.particleNumber(T_gas)
ofile = open("gas_{0}.dat".format(sys.argv[2][-4:]), "wb")
nin = 0
for gaspart, j in zip(data.iterator(T_gas), range(ngas)):
    if gaspart.x() > -100. and gaspart.x() < 100.:
      if gaspart.y() > -100. and gaspart.y() < 100.:
        if gaspart.z() > -100. and gaspart.z() < 100.:
#          ofile.write(str(gaspart.x()) + "\t" + str(gaspart.y()) + "\t" + str(gaspart.z()) + "\t" + str(gaspart.density()) + "\t" + str(gaspart.sml()) + "\n")
          ofile.write(struct.pack('f', gaspart.x()))
          ofile.write(struct.pack('f', gaspart.y()))
          ofile.write(struct.pack('f', gaspart.z()))
          ofile.write(struct.pack('f', gaspart.density()))
          ofile.write(struct.pack('f', gaspart.sml()))
          ofile.write(struct.pack('f', gaspart.temperature()))
          nin += 1
ofile.close()

print "extracted", nin, "from", ngas, "gas particles"

ofile = open("stars_{0}.dat".format(sys.argv[2][-4:]), "wb")
nstar = data.particleNumber(T_star)
for starpart, j in zip(data.iterator(T_star), range(nstar)):
    ofile.write(struct.pack('f', starpart.x()))
    ofile.write(struct.pack('f', starpart.y()))
    ofile.write(struct.pack('f', starpart.z()))
ofile.close()
