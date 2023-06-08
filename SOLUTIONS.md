# CASA VLBI calibration tutorial

- Prepare your data by appending the data necessary for the a-priori
  (amplitude) calibration.

  `$ casa --no-gui -c ~/src/casa-vlbi/append_tsys.py n14c3.antab n14c3_1_1.IDI?`

  `$ casa --no-gui -c ~/src/casa-vlbi/append_gc.py n14c3.antab n14c3_1_1.IDI1`

- Load the data into CASA (converting it into MeasurementSet format):

  `>>> importfitsidi(vis='n14c3.ms', fitsidifile=['n14c3_1_1.IDI1', 'n14c3_1_1.IDI2'], constobsid=True, scanreindexgap_s=15)`

- Use the listobs task to get an overview of your observation

  Hint: you can use the "listfile" parameter of the listobs task to
  write the output to a file for future reference.

  `>>> listobs(vid='n14c3.ms')`
  `>>> listobs(vis='n14c3.ms', listfile='n14c3.listobs')`

- Check that the Tsys measurements and the gain curves made it into
  your data MeasurementSet.

  Hint: The Tsys measurements are stored in the SYSCAL sub-table of your
  MeasurementSet.  The gain-curves are stored in the GAIN_CURVE
  sub-table.  You can use the CASA table tool to take a detailed look at
  these tables and check for example the number of rows in the table (to
  check that the table isn't empty) or write the contents of the table
  to an (ASCII) text file.

  `>>> tb.open('n14c3.ms/SYSCAL')`
  `>>> tb.nrows()`
  `40440`
  `>>> tb.toasciifmt('SYSCAL.txt')`

  `>>> tb.open('n14c3.ms/GAIN_CURVE')`
  `>>> tb.nrows()`
  96
  `>>> tb.toasciifmt('GAIN_CURVE.txt')`

- Save the current flagging state.

  Hint: This can be done using the flagmanager task.

  `>>> flagmanager(vis='n14c3.ms', mode='save', versionname='initial')`

- Plot the amplitudes of your data with the plotms task to get some idea
  about what's there.  For example, plot amplitude as a function of
  time.

  Hint: You probably don't want to plot all the data in the data set.
  You can use the "scan" parameter to plot data for only single scan.

  `>>> plotms(vis="n14c3.ms", xaxis="time",yaxis="amp", scan='2', coloraxis='baseline')`


- Use the gencal task to generate an amplitude calibration table based
  on the Tsys measurements, but don't specify the uniform parameter.

  `>>> gencal(vis='n14c3.ms', caltable='n14c3.tsys_uniform', caltype='tsys')`

- Apply the calibration to your data.

  `>>> applycal(vis='n14c3.ms', gaintable='n14c3.tsys_uniform')`

- Plot your data again using the same command as you used before.
  What differences do you see?

  `>>> plotms(vis="n14c3.ms", xaxis="time", yaxis="amp", scan='2', coloraxis='baseline')`

- Look at the CASA log file.  There should be a line like this:

  2023-05-30 12:35:34     INFO    applycal::::casa           B TSYS: In: 27607552 / 236298240   (11.683350667360028%) --> Out: 234054528 / 236298240   (99.05047451898076%) (n14c3.tsys)

- Can you explain what happened to your data?

- Since applying the calibration changed the flags, restore the flags
  that you saved before.

  Hint: This can also be done using the flagmanager task.

  `>>> flagmanager(vis='n14c3.ms', mode='restore', versionname='initial')`

- Now regenerate the Tsys calibration table but add the "unifrom=False" parameter.

  `>>> gencal(vis='n14c3.ms', caltable='n14c3.tsys', caltype='tsys', uniform=False)`

- Apply the calibration to your data.

  `>>> applycal(vis='n14c3.ms', gaintable='n14c3.tsys')`

- Plot your data once more using the same command as before.  Do
  things look better now?

  `>>> plotms(vis="n14c3.ms", xaxis="time", yaxis="amp", scan='2', coloraxis='baseline')`

- You may want to look at the CASA log file again.

  2023-05-30 13:28:52     INFO    applycal::::casa           B TSYS: In: 27607552 / 236298240   (11.683350667360028%) --> Out: 27951968 / 236298240   (11.829105455884902%) (n14c3.tsys2)

- Use the gencal task to generate a calibration table for the gain curves and DPFU.

  `>>> gencal(vis='n14c3.ms', caltable='n14c3.gcal', caltype='gc')`

- Apply both the Tsys and gain curve calibration.

  `>>> applycal(vis='n14c3.ms', gaintable=['n14c3.tsys', 'n14c3.gcal'])`

- Now make a more detailed amplitude vs. time plot of your data.  Pick a
  reference station (e.g. EF) using the "refant" parameter, pick one of
  the parallel polarizations (e.g. LL) using the "correlation"
  parameter, pick a spectral window (e.g. 2) using the "spw" parameter
  and average the different channels for more signal-to-noise using the
  "avgchannel" paremeter.  Consult the listobs output for the number of
  channels to average.

  First plot the uncalibrated data, and then plot the calibrated
  data.  This can be done using the "ydatacolumn" paraeter.

  `>>> plotms(vis="n14c3.ms", xaxis="time",yaxis="amp", scan='2', coloraxis='baseline', ydatacolumn='data', correlation='ll', avgchannel="32", scalar=True, spw='2',antenna='EF')`

  `>>> plotms(vis="n14c3.ms", xaxis="time",yaxis="amp", scan='2', coloraxis='baseline', ydatacolumn='corrected', correlation='ll', avgchannel="32", scalar=True, spw='2',antenna='EF')`


- Now plot the weights instead of the amplitude.  What has happened to the weights  during calibration?

  `>>> plotms(vis="n14c3.ms", xaxis="time",yaxis="wt", scan='2', coloraxis='baseline', ydatacolumn='data', correlation='ll', avgchannel="32", scalar=True, spw='2',antenna='EF')`

  `>>> plotms(vis="n14c3.ms", xaxis="time",yaxis="wt", scan='2', coloraxis='baseline', ydatacolumn='corrected', correlation='ll', avgchannel="32", scalar=True, spw='2',antenna='EF')`


Before calibration the weights were (nearly) identical for all the
baselines.  After calibration the weight are different.  The weight of
baselines to less sensitive antennas will be lower than those to more
sesitive antennas
