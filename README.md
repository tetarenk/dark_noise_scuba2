# dark_noise_scuba2
Python script to characterize bolometer noise in the SCUBA-2 instrument on the JCMT.

## General Description
The script grabs the effective NEPs for all SCUBA-2 sub-arrays from the JCMT logs, and writes them into output files for 450/850 um. Time series plots of effective NEPs vs time, and effective NEP histograms separated by sub-array are also output.

## Requires
* `numpy`
* `matplotlib`
* `astropy`

## Notes
* Needs to be run on an EAO computer.
* Compatible with Python 3.
