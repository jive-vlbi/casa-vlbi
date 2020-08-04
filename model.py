#! /usr/bin/env python
#
# Copyright (C) 2018
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
import json
import optparse
import os
import struct
import sys
import time as tm
import urlparse

import pyfits
import numpy as np
import scipy as sp
from scipy import fftpack, signal

from vex import Vex

def vex2time(str):
    tupletime = tm.strptime(str, "%Yy%jd%Hh%Mm%Ss");
    return tm.mktime(tupletime)

class ModelTable:
    def __init__(self):
        self.time = []
        self.time_interval = []
        self.source_id = []
        self.antenna_no = []
        self.array = []
        self.freqid = []
        self.gdelay = []
        return
    pass

class ScanInfo:
    def __init__(self, vex, station, scan):
        self.scan = scan

        source = vex['SCHED'][scan]['source']
        self.source = vex['SOURCE'][source]['source_name']
        self.start = vex2time(vex['SCHED'][scan]['start'])
        self.length = 0
        for transfer in vex['SCHED'][scan].getall('station'):
            if transfer[0] == station:
                self.length = int(transfer[2].split()[0])
                pass
            continue

        clock = vex['STATION'][station]['CLOCK']
        clock_early = vex['CLOCK'][clock]['clock_early']
        self.clock = np.zeros(2)
        self.clock_epoch = vex2time(clock_early[2])
        value = clock_early[1].split()
        self.clock[0] = float(value[0])
        if value[1] == 'usec':
            self.clock[0] *= 1e-6
            pass
        value = clock_early[3].split()
        # If clock rate unit is missing, assume usec/sec.
        if len(value) == 1:
            value.append('usec/sec')
            pass
        self.clock[1] = float(value[0])
        if value[1] == 'usec/sec':
            self.clock[1] *= 1e-6
            pass
        return
    pass

delay_header_v0 = "=I2sx"
delay_header_v1 = "=II2sx"
delay_scan = "=80sx"
delay_source = "=80sxI"
delay_entry = "=7d"

def parse_model(info, delay_file):
    fp = open(delay_file, 'r')
    buf = fp.read(struct.calcsize(delay_header_v0))
    hdr = struct.unpack(delay_header_v0, buf)
    if hdr[0] > 3:
        fp.seek(0)
        buf = fp.read(struct.calcsize(delay_header_v1))
        hdr = struct.unpack(delay_header_v1, buf)
        version = hdr[1]
    else:
        version = -1
        pass

    scan = info.scan

    if version == 0:
        buf = fp.read(struct.calcsize(delay_scan))
        scan = struct.unpack(delay_scan, buf)[0].strip('\0')
        pass
    buf = fp.read(struct.calcsize(delay_source))
    src = struct.unpack(delay_source, buf)
    source = src[0].strip()
    mjd = src[1]
    start = None
    t = []
    d = []
    u = []; v = [];  w = []
    while fp:
        buf = fp.read(struct.calcsize(delay_entry))
        if not buf:
            break
        if len(buf) < struct.calcsize(delay_entry):
            break
        delay = struct.unpack(delay_entry, buf)
        if delay[0] == 0 and delay[4] == 0:
            if (source == info.source and
                scan == info.scan and
                start >= info.start and
                start < info.start + info.length):
                return (t, d, u, v, w)
            if version == 0:
                buf = fp.read(struct.calcsize(delay_scan))
                if not buf:
                    break
                scan = struct.unpack(delay_scan, buf)[0].strip('\0')
                pass
            buf = fp.read(struct.calcsize(delay_source))
            if not buf:
                break
            src = struct.unpack(delay_source, buf)
            source = src[0].strip()
            mjd = src[1]
            start = None
            t = []
            d = []
            u = []; v = [];  w = []
            continue
        if start == None:
            start = (mjd - 40587) * 86400 + delay[0]
            pass
        t.append(delay[0])
        u.append(delay[1])
        v.append(delay[2])
        w.append(delay[3])
        d.append(delay[4])
        continue
    return

def create_splines(interval, x, y):
    times = []
    splines = []
    diff_max = 0.0
    x = np.array(x)
    y = np.array(y)
    akima = sp.interpolate.Akima1DInterpolator(x, y)
    while len(x) > 1:
        points = min(interval + 1, len(x))
        u = np.linspace(0, points - 1, 7)
        v = np.array(map(lambda x: akima(x), x[0] + u))
        z = np.polyfit(u, v, 5)
        times.append(x[0])
        splines.append(np.flip(z, 0))
        poly = np.poly1d(z)
        a = np.arange(0, points - 1, 0.125)
        b = np.array(map(lambda x: poly(x), a))
        q = np.arange(x[0], x[0] + points - 1, 0.125)
        r = np.array(map(lambda x: akima(x), q))
	diff = np.max(np.absolute(r - b))
        if diff > diff_max:
            diff_max = diff
            pass
        x = x[points - 1:]
        y = y[points - 1:]
        continue
    return (times, splines, diff_max)

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

class IdiData:
    def __init__(self, idifiles):
        hdulist = pyfits.open(idifiles[0])

        # Fundamentals
        header = hdulist['ARRAY_GEOMETRY'].header
        self.obscode = header['OBSCODE']
        self.n_stokes = header['NO_STKD']
        self.stk_1 = header['STK_1']
        self.n_band = header['NO_BAND']
        self.n_chan = header['NO_CHAN']
        self.ref_freq = header['REF_FREQ']
        self.chan_bw = header['CHAN_BW']
        self.ref_pixl = header['REF_PIXL']
        self.rdate = header['RDATE']

        # Create antenna mapping
        tbhdu = hdulist['ARRAY_GEOMETRY']
        assert tbhdu.header['EXTVER'] == 1
        antennas = [s.strip() for s in tbhdu.data['ANNAME']]
        self.antenna_map = dict(zip(antennas, tbhdu.data['NOSTA']))

        # Create list of frequencies
        tbhdu = hdulist['FREQUENCY']
        assert tbhdu.data['FREQID'][0] == 1
        freqs = tbhdu.data['BANDFREQ'][0]
        self.freqs = np.array([x + self.ref_freq for x in freqs])
        assert len(self.freqs) == self.n_band

        # Create source mapping
        tbhdu = hdulist['SOURCE']
        if 'ID_NO.' in tbhdu.columns.names:
            source_id_col = 'ID_NO.'
        else:
            source_id_col = 'SOURCE_ID'
            pass
        self.source_map = dict(zip(tbhdu.data['SOURCE'],
                                   tbhdu.data[source_id_col]))

        # Reference date in convenient form
        tupletime = tm.strptime(self.rdate, "%Y-%m-%d")
        self.reftime = tm.mktime(tupletime)

        # Create an index
        self.idx = []
        source_id = -1
        self.last_time = float("-inf")
        self.first_time = float("inf")
        for idifile in idifiles:
            hdulist = pyfits.open(idifile)
            tbhdu = hdulist['UV_DATA']
            if 'SOURCE' in tbhdu.columns.names:
                source_id_col = 'SOURCE'
            else:
                source_id_col = 'SOURCE_ID'
                pass
            for data in tbhdu.data:
                jd = data['DATE']
                time = (jd - 2440587.5 + data['TIME']) * 86400
                if time > self.last_time:
                    self.last_time = time
                    pass
                if time < self.first_time:
                    self.first_time = time
                    pass
                if data[source_id_col] != source_id:
                    source_id = data[source_id_col]
                    self.idx.append((time, source_id))
                    pass
                continue
            continue
        self.idx.sort()
        return

    def find_source(self, time):
        if time < self.first_time or time > self.last_time:
            return -1
        source_id = -1
        i = 0
        try:
            while time >= self.idx[i][0]:
                source_id = self.idx[i][1]
                i += 1
                continue
        except:
            pass
        return source_id


os.environ['TZ'] = 'UTC'
tm.tzset()

usage = "usage %prog [options] vexfile ctrfile fitsfile"
parser = optparse.OptionParser(usage=usage)
(options, args) = parser.parse_args(sys.argv[1:])

def process_job(vex, ctrl, idi, data, interval=120):
    exper_name = ctrl['exper_name']

    delay_uri = ctrl['delay_directory']
    delay_directory = urlparse.urlparse(delay_uri).path
    delay_directory = delay_directory[17:]

    scan = ctrl['scans'][0]

    for station in ctrl['stations']:
        info = ScanInfo(vex, station, scan)
        antenna = station.upper()

        delay_file = exper_name + '_' + station + '.del'
        delay_file = os.path.join(delay_directory, delay_file)
        (t, d, u, v, w) = parse_model(info, delay_file)
        (t, d, diff_max) = create_splines(interval, t, d)
        data.time_interval.extend(len(t) * [float(interval) / 86400])
        data.source_id.extend(len(t) * [idi.source_map[info.source]])
        data.antenna_no.extend(len(t) * [idi.antenna_map[antenna]])
        data.array.extend(len(t) * [1])
        data.freqid.extend(len(t) * [1])
        start = info.start
        for spline in d:
            data.time.append(float(start - idi.reftime) / 86400)
            spline[0] += info.clock[0] + (start - info.clock_epoch) * info.clock[1]
            spline[1] += info.clock[1]
            start += interval
            continue
        data.gdelay.extend(d)
        continue
    return

lisfile = args[0]
idifiles = args[1:]

fp = open(lisfile, 'r')
line = fp.readline()
vexfile = line.split()[1]

ctrlfiles = []
for line in fp:
    if line.startswith('+'):
        corfile = line.split()[1]
        ctrlfile = os.path.splitext(corfile)[0] + ".ctrl"
        ctrlfiles.append(ctrlfile)
        pass
    continue

vex = Vex(vexfile)

idi = IdiData(idifiles)
pols = ['R']
if idi.n_stokes > 1:
    pols.append('L')
    pass

data = ModelTable()

for ctrlfile in ctrlfiles:
    fp = open(ctrlfile, 'r')
    ctrl = json.load(fp)
    fp.close()
    process_job(vex, ctrl, idi, data)
    continue

time = np.array(data.time)
time_interval = np.array(data.time_interval)

source_id = np.array(data.source_id)
antenna_no = np.array(data.antenna_no)
array = np.array(data.array)
freqid = np.array(data.freqid)

i_far_rot = np.zeros(len(time))
freq_var = np.zeros(len(time))

n_poly = 6
nrows = len(data.gdelay)
freqs = np.repeat(idi.freqs, n_poly)

gdelay = np.array(data.gdelay)
gdelay = np.repeat(gdelay, idi.n_band, axis=0)
gdelay = gdelay.reshape(nrows, n_poly * idi.n_band)

pdelay = map(lambda x: x * freqs, gdelay)

grate = []
for poly in data.gdelay:
    r = np.zeros(n_poly)
    for n in range(n_poly - 1):
        r[n] = ((n + 1) * poly[n + 1])
        continue
    grate.append(r)
grate = np.array(grate)
grate = np.repeat(grate, idi.n_band, axis=0)
grate = grate.reshape(nrows, n_poly * idi.n_band)

prate = map(lambda x: x * freqs, grate)

disp = np.zeros(len(time))
ddisp = np.zeros(len(time))

hdulist = pyfits.open(idifiles[0])
header = hdulist['ARRAY_GEOMETRY'].header
obscode = header['OBSCODE']
n_stokes = header['NO_STKD']
stk_1 = header['STK_1']
n_band = header['NO_BAND']
n_chan = header['NO_CHAN']
ref_freq = header['REF_FREQ']
chan_bw = header['CHAN_BW']
ref_pixl = header['REF_PIXL']
rdate = header['RDATE']

tbhdu = hdulist['ARRAY_GEOMETRY']
assert tbhdu.header['EXTVER'] == 1
antennas = [s.strip() for s in tbhdu.data['ANNAME']]
antenna_map = dict(zip(antennas, tbhdu.data['NOSTA']))

tbhdu = hdulist['FREQUENCY']
assert tbhdu.data['FREQID'][0] == 1
freqs = tbhdu.data['BANDFREQ'][0]
freqs = np.array(map (lambda x: x + ref_freq, freqs))
assert len(freqs) == n_band

tbhdu = hdulist['SOURCE']
if 'ID_NO.' in tbhdu.columns.names:
    source_id_col = 'ID_NO.'
else:
    source_id_col = 'SOURCE_ID'
    pass
source_map = dict(zip(tbhdu.data['SOURCE'], tbhdu.data[source_id_col]))

tupletime= tm.strptime(rdate, "%Y-%m-%d")
ref = tm.mktime(tupletime)

cols = []
col = pyfits.Column(name='TIME', format='1D', unit='DAYS', array=time)
cols.append(col)
col = pyfits.Column(name='TIME_INTERVAL', format='1E', unit='DAYS', array=time_interval)
cols.append(col)
col = pyfits.Column(name='SOURCE_ID', format='1J', array=source_id)
cols.append(col)
col = pyfits.Column(name='ANTENNA_NO', format='1J', array=antenna_no)
cols.append(col)
col = pyfits.Column(name='ARRAY', format='1J', array=array)
cols.append(col)
col = pyfits.Column(name='FREQID', format='1J', array=freqid)
cols.append(col)
col = pyfits.Column(name='I.FAR.ROT', format='1E', unit='RAD/METER**2', array=i_far_rot)
cols.append(col)
format = '%dE' % n_band
col = pyfits.Column(name='FREQ.VAR', format=format, unit='HZ', array=freq_var)
cols.append(col)
format = '%dD' % (n_poly * n_band)
col = pyfits.Column(name='PDELAY_1', format=format, unit='TURNS', array=pdelay)
cols.append(col)
col = pyfits.Column(name='GDELAY_1', format=format, unit='SECONDS', array=gdelay)
cols.append(col)
col = pyfits.Column(name='PRATE_1', format=format, unit='HZ', array=pdelay)
cols.append(col)
col = pyfits.Column(name='GRATE_1', format=format, unit='SEC/SEC', array=grate)
cols.append(col)
col = pyfits.Column(name='DISP_1', format='1E', unit='SECONDS', array=disp)
cols.append(col)
col = pyfits.Column(name='DDISP_1', format='1E', unit='SEC/SEC', array=ddisp)
cols.append(col)
if len(pols) > 1:
    col = pyfits.Column(name='PDELAY_2', format=format, unit='TURNS', array=pdelay)
    cols.append(col)
    col = pyfits.Column(name='GDELAY_2', format=format, unit='SECONDS', array=gdelay)
    cols.append(col)
    col = pyfits.Column(name='PRATE_2', format=format, unit='HZ', array=prate)
    cols.append(col)
    col = pyfits.Column(name='GRATE_2', format=format, unit='SEC/SEC', array=grate)
    cols.append(col)
    col = pyfits.Column(name='DISP_2', format='1E', unit='SECONDS', array=disp)
    cols.append(col)
    col = pyfits.Column(name='DDISP_2', format='1E', unit='SEC/SEC', array=ddisp)
    cols.append(col)
    pass
coldefs = pyfits.ColDefs(cols)

header = pyfits.Header()
header['EXTNAME'] = 'INTERFEROMETER_MODEL'
header['EXTVER'] = 1
header['TABREV'] = 2
header['OBSCODE'] = idi.obscode
header['NO_STKD'] = idi.n_stokes
header['STK_1'] = idi.stk_1
header['NO_BAND'] = idi.n_band
header['NO_CHAN'] = idi.n_chan
header['REF_FREQ'] = idi.ref_freq
header['CHAN_BW'] = idi.chan_bw
header['REF_PIXL'] = ref_pixl
# Repeat the reference data even though the FITS-IDI standard doesn't
# seem to require it.
header['RDATE'] = idi.rdate
header['NPOLY'] = n_poly
header['NO_POL'] = len(pols)

tbhdu = pyfits.BinTableHDU.from_columns(coldefs, header)

hdulist = pyfits.open(idifiles[0], mode='append')
hdulist.append(tbhdu)
hdulist.close()
