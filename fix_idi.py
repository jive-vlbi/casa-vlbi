#! /usr/bin/env python
#
# Copyright (C) 2020
# Joint Institute for VLBI ERIC, Dwingeloo, The Netherlands
#
# This library is free software; you can redistribute it and/or modify it
# under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
# License for more details.
#
# You should have received a copy of the GNU Library General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
#
import optparse

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

class IdiHDU(pyfits.PrimaryHDU):
    @classmethod
    def match_header(cls, header):
        try:
            keyword = header.cards[0].keyword
        except:
            keyword = header.ascard[0].key
            pass
        return (keyword == 'SIMPLE' and 'GROUPS' in header and
                header['GROUPS'] == True and 'NAXIS' in header and
                header['NAXIS'] == 0)

pyfits.register_hdu(IdiHDU)

usage = "usage %prog [options] infile outfile"
parser = optparse.OptionParser(usage=usage)
(options, args) = parser.parse_args()
if len(args) != 2:
    parser.error("incorrect number of arguments")

hdulist = pyfits.open(args[0])
header = hdulist['UV_DATA'].header
tbhdu = hdulist['UV_DATA']

data = tbhdu.data
mask = data['baseline'] != 0
newdata = data[mask]
newhdu = pyfits.BinTableHDU(data = newdata, header = header, name ='UV_DATA')

hdus = []
for hdu in hdulist:
    if hdu == hdulist['UV_DATA']:
        continue
    hdus.append(hdu)
    continue
hdus.append(newhdu)

newhdulist = pyfits.HDUList(hdus)
newhdulist.writeto(args[1])
newhdulist.close()

hdulist.close()
