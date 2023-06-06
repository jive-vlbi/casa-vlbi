# CASA VLBI calibration tutorial

- Prepare your data by appending the data necessary for the a-priori
  (amplitude) calibration.

- Load the data into CASA (converting it into MeasurementSet format):

- Use the listobs task to get an overview of your observation

  Hint: you can use the "listfile" parameter of the listobs task to
  write the output to a file for future reference.

- Check that the Tsys measurements and the gain curves made it into your
  data MeasurementSet.

  Hint: The Tsys measurements are stored in the SYSCAL sub-table of your
  MeasurementSet.  The gain-curves are stored in the GAIN_CURVE
  sub-table.  You can use the CASA table tool to take a detailed look at
  these tables and check for example the number of rows in the table (to
  check that the table isn't empty) or write the contents of the table
  to an (ASCII) text file.

- Save the current flagging state.

  Hint: This can be done using the flagmanager task.

- Plot the amplitudes og your data with the plotms task to get some idea
  about what's there.  For example, plot amplitude as a function of
  time.

  Hint: You probably don't want to plot all the data in the data set.
  You can use the "scan" parameter to plot data for only single scan.

- Use the gencal task to generate an amplitude calibration table based
  on the Tsys measurements, but don't specify the uniform parameter.

- Apply the calibration to your data.

- Plot your data again using the same command as you used before.
  What differences do you see?

- Look at the CASA log file.  There should be a line like this:

  2023-05-30 12:35:34     INFO    applycal::::casa           B TSYS: In: 27607552 / 236298240   (11.683350667360028%) --> Out: 234054528 / 236298240   (99.05047451898076%) (n14c3.tsys)

  Can you explain what happened to your data?

- Since applying the calibration changed the flags, restore the flags
  that you saved before.

  Hint: This can also be done using the flagmanager task.

- Now regenerate the Tsys calibration table but add the
  "unifrom=False" parameter.

- Apply the calibration to your data.

  Hint: this can be done using the applycal task

- Plot your data once more using the same command as before.  Do
  things look better now?

  You may want to look at the CASA log file again.

  2023-05-30 13:28:52     INFO    applycal::::casa           B TSYS: In: 27607552 / 236298240   (11.683350667360028%) --> Out: 27951968 / 236298240   (11.829105455884902%) (n14c3.tsys2)

- Use the gencal task to generate a calibration table for the gain
  curves and DPFU.

- Apply both the Tsys and gain curve calibration.

- Now make a more detailed amplitude vs. time plot of your data.  Pick a
  reference station (e.g. EF) using the "refant" parameter, pick one of
  the parallel polarizations (e.g. LL) using the "correlation"
  parameter, pick a spectral window (e.g. 2) using the "spw" parameter
  and average the different channels for more signal-to-noise using the
  "avgchannel" paremeter.  Consult the listobs output for the number of
  channels to average.

- First plot the uncalibrated data, and then plot the calibrated
  data.  This can be done using the "ydatacolumn" paraeter.

- Now plot the weights instead of the amplitude.  What has happened to
  the weights during calibration?

  Hint: the weights can be plotted by using yaxis="wt" in the plotms
  command.
