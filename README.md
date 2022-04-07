# CASA-VLBI Scripts

Scripts to assist VLBI data reduction in [CASA](https://casa.nrao.edu).
These scripts focus on the a-priori (initial) gain calibration of the
different VLBI antennas, and the a-priori flagging commands, both provided
by the EVN Data Archive.

All scripts can run as a standalone program or being called within the
CASA environment.


## A-priori gain calibration

The a-priori gain calibration sets the initial gains of each antenna
based on their gain curves and the recorded system temperatures during
the observation.

To prepare the data in order to be able to apply the gain calibration in
CASA you will need to run two steps before the FITS-IDI files are converted
into a Measuarement Set (MS) format.


### Appending the system temperatures to the FITS IDI files

You can run it either via the standalone script:

```bash
append_tsys.py  {antabfile}  {idifiles}
```
or via the Python module:
```python
casavlbitools.fitsidi.append_tsys(antabfile, idifiles)
```

This will append Tsys measurements from the file named by `antabfile` (which
should be in AIPS ANTAB format or VLBA VLOG format) to the FITS-IDI
file named by `idifiles`.  If an observation is split into multiple
FITS-IDI files, `idifile` should be a list of file names.



The functionality of the append_tsys.py, gc.py and flag.py scripts is
also available in Python modules.  The following interfaces are
available:


### Appending the Gain Curve and sensitivity to the FITS IDI files

You can run it either via the standalone script:
```bash
append_gc.py  {antabfile}  {idifiles}
```
or via the Python module:
```python
casavlbitools.fitsidi.append_gc(antabfile, idifile)
```
This will append gain curves and sensitivity information from the file named by
`antabfile` (which should be in AIPS ANTAB format) to the FITS-IDI
file named by `idifile`.


### Create the Gain Curve calibration table

*-- Note that if you already appended the system temperatures and the gain curves into the
FITS-IDI files before convert it to a MS file, you can directly run `gencal` in CASA
with such information without requiring this step. --*

To create an external gain curve calibration table with the responses for each antenna:
```bash
gc.py {antabfile}  {gc}  --min-elevation=0.0  --max=elevation=90.0
```
or via the Python module:
```python
casavlbitools.casa.convert_gaincurve(antabfile, gc, min_elevation=0.0, max_elevation=90.0)
```
This creates a gaincurve table with name `gc` from gain curves from the file
named by `antablfile` (which should be in AIPS ANTAB format or VLBA
gains format).  The gain curves are sampled between `min_elevation`
and `max_elevation`, which should be specified in degrees, and then
refitted.


## Convert the file format or the a-priori flagging

If you have a flag file in AIPS UVFLG format, you can convert it into a
CASA FLAG format with the following script:
```bash
flag.py  {infile}  {idifiles}
```
or using the Python module:
```python
casavlbitools.fitsidi.convert_flags(infile, idifiles, outfp=sys.stdout, outfile=None)
```

This converts the flag file named by `infile` (which should be in AIPS UVFLG
format) into a format that can be used by CASA.  The name of the
FITS-IDI files corresponding to the observation should be passed in
`idifiles`.  By default the converted output is written to stdout.
Output is written to a file named `outpufile` instead if that
parameter is provided.
