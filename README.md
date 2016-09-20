# OrbitalDetermination35

Orbital Determination via the Method of Gauss

Alexander Davenport

Original version 07/22/2016, SSP NMT 2016

Current version: 2.0.0 (09/19/2016)

##Usage

Run the code from the command line like so: 

`python main.py <input-file-path>`

You will be prompted to pick a root by the program after a second. Pick a distance (in AU) that is reasonable for the class of asteroid you are observing. Wait for the code to run through the iteration and print out your data at the end if you feel the need. 

##Input file format

The format of the input file is very specific in order but is tolerant in some regards, like the number of decimal places in the sun-earth vector components, time, or coordinates.

An example input file might look like:
```
10150 (1994 PN)
2016-06-22 08:03:10.429	17:26:58.87 30:46:11.30 -0.0195340791840548 0.9323384266123208 0.4041408400938959
2016-06-30 07:19:24.379	17:08:12.58 27:19:02.51 -0.1537983334829756 0.9221125109705270 0.3997035045290656
2016-07-07 05:21:27.350 16:51:52.97 23:04:41.48 -0.2682470588825682 0.8998279354659114 0.3900449087482382
2016-07-13 03:42:48.670 16:38:37.64 18:31:04.40 -0.3635214279037087 0.8710218039698151 0.3775632429881836
```
With the format of each line being:
```
Date       Time       RA(HMS)     Dec(DMS)    Sun-Earth Vector (AU)
YYYY-MM-DD HH:MM:SS.s HH:MM:SS.ss DD:MM:SS.ss X Y Z
```

##Further Reading
For more information about the Method of Gauss and the calculations contained within, please consult [this document](https://www.overleaf.com/read/ywydpqpzxhnk).
