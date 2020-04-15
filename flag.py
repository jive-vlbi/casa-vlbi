#! /usr/bin/env python
#
# Copyright (C) 2016, 2017
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
import datetime
import inspect
import optparse
import os
import sys

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

filename = inspect.getframeinfo(inspect.currentframe()).filename
sys.path.append(os.path.dirname(os.path.realpath(filename)))
from casavlbitools import key

def find_first_dobs(idifiles):
    first_dobs = datetime.datetime(datetime.MAXYEAR, 12, 31)
    for match in idifiles:
        hdulist = pyfits.open(match)
        dobs = hdulist['UV_DATA'].header['DATE-OBS']
        dobs = datetime.datetime.strptime(dobs, '%Y-%m-%d')
        if dobs < first_dobs:
            first_dobs = dobs
    return first_dobs

def convert_flags(infile, idifiles, outfp=sys.stdout, outfile=None):
    # Make sure we're dealing with a list
    if not isinstance(idifiles, list):
        idifiles = [idifiles]
        pass

    # Create new output file if requested to do so
    if outfile:
        outfp = open(outfile, 'w')
        pass

    first_date = find_first_dobs(idifiles)

    infp = open(infile, 'r')
    flags = key.read_keyfile(infp)
    infp.close()

    flags = [dict(x) for x in flags]

    for flag in flags:
        antenna = flag['ANT_NAME'].upper()
        timerange = flag['TIMERANG']
        if timerange == [0.0, 0.0, 0.0, 0.0, 400.0, 0.0, 0.0, 0.0]:
            timerange = ""
        else:
            year = datetime.datetime(first_date.year, 1, 1)
            date1 = year + datetime.timedelta(timerange[0] - 1)
            date1 = date1.strftime("%Y/%m/%d")
            date2 = year + datetime.timedelta(timerange[4] - 1)
            date2 = date2.strftime("%Y/%m/%d")
            timerange = "%s/%02d:%02d:%02d~%s/%02d:%02d:%02d" % \
                (date1, timerange[1], timerange[2], timerange[3],
                 date2, timerange[5], timerange[6], timerange[7])
            pass
        if 'BIF' in flag:
            spw = "%d~%d" % (flag['BIF'] -1, flag['EIF'] - 1)
        else:
            spw = ""
            pass
        if 'BCHAN' in flag:
            spw = spw + ":%d~%d" % (flag['BCHAN']-1, flag['ECHAN']-1)
            pass
        cmd = "antenna='%s'" % antenna
        if timerange:
            cmd = cmd + " timerange='%s'" % timerange
            pass
        if spw:
            cmd = cmd + " spw='%s'" % spw
            pass
        outfp.write(cmd + "\n")
    return

if __name__ == "__main__":
    # Allow this scrupt to be invoked using casa -c
    try:
        i = sys.argv.index("-c") + 2
    except:
        i = 1
        pass

    usage = "usage %prog [options] uvflagfile fitsfile"
    parser = optparse.OptionParser(usage=usage)
    (options, args) = parser.parse_args(sys.argv[i:])
    convert_flags(args[0], args[1:])
