# MIRI-LRS-dithers

Code to generate and visualize the LRS commissioning &amp; calibration dither patterns

## Contents of this repository

This repo contains code and documentation to generate dither patterns for the JWST-MIRI Low Resolution Spectrometer (LRS). The majority of these patterns are designed to support the commissioning and calibration activities, i.e. not for general science use. 

The output of the code is 1 or more dither pattern tables that can be delivered to the Astronomers Proposal Tool (APT) team for ingestion into APT and into the Project Reference Database (PRD). 


## Dependencies

* numpy
* matplotlib (for visualization of the patterns)
* astropy
* pySIAF (part of astroconda)


## Key to dither patterns

The tables below matches dither patterns to commissioning activity requests (CARs) and calibration activities (CALs).

### 1. Slit


| Pattern name          | No. positions    | CAR      | CAL      |
|-----------------------|------------------|----------|----------|
|LRS-1pixscan-fullslit  |         45       |          |          |
|LRS-2pixscan-fullslit  |         23       |          |          |
|LRS-7pixscan-fullslit  |          7       |          |          |
|LRS-7x3-center         |         21       |          |          |
|LRS-7x3-nod1           |         21       |          |          |
|LRS-7x3-nod2           |         21       |          |          |
|LRS-intrapixel-center  |         11       |          |          |
|LRS-intrapixel-nod1    |         11       |          |          |
|LRS-intrapixel-nod2    |         11       |          |          |
|LRS-long-cross-center  |         17       |          |          |
|LRS-long-cross-nod1    |         17       |          |          |
|LRS-long-cross-nod2    |         17       |          |          |
|LRS-short-cross-center |          9       |          |          |
|LRS-short-cross-nod1   |          9       |          |          |
|LRS-short-cross-nod2   |          9       |          |          |


### 2. Slitless


| Pattern name          | No. positions    | CAR      | CAL      |
|-----------------------|------------------|----------|----------|
|LRS-1pix-slitless-long |       67         |          |          |
|LRS-1pix-slitless-short|       51         |          |          |
|LRS-2pix-slitless-short|       27         |          |          |
|LRS-7x3-center-slitless|       21         |          |          |
|LRS-7pix-9x3-slitless  |       27         |          |          |
|LRS-5pix-8x4-slitless  |       32         |          |          |
|LRS-intrapixel-center  |       11         |          |          |
|LRS-7pixscan-slitless  |       9          |          |          |
|LRS-5pixscan-slitless  |       8          |          |          |

