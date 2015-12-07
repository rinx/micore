MICO-RE
===
MICO-RE 【miko-li】  
Minimal Implementation of Cloud Optical Retrieval

MICO-RE estimates cloud optical thickness (COT) and
cloud droplet effective radius (CDER) from reflectances in two wavelengths
by using look-up table (Nakajima and King, 1990) and Gauss-Newton method (Rodgers, 2000).  
Akima method (Akima, 1970) is used as a method of interpolation.


Cloud Optical Retrieval Code
---

### Installation

Just execute `make` command in `src` directory.

    $ cd src
    $ make

You can change your compiler setting in the `Makefile` (The default compiler is `gfortran`).

### Usage

    $ ./micore [lutfile] [reflectance1] [reflectance2]

like following:

    $ ./micore ../example/lut_860_2130.bin 0.553 0.343


Look-up table
---
Look-up table used in this program should be like following:

ex.)

| COT | CDER | REF1 | REF2 |
|:----|:-----|:-----|:-----|
| 1.0 | 3.0  | 0.20 | 0.10 |
| 1.0 | 5.0  | 0.20 | 0.13 |
| 1.0 | 7.0  | 0.21 | 0.15 |
| 1.0 | 9.0  | 0.23 | 0.14 |
| ... | ...  | ...  | ...  |
| ... | ...  | ...  | ...  |
| 2.0 | 3.0  | 0.26 | 0.14 |
| 2.0 | 5.0  | 0.26 | 0.18 |
| 2.0 | 7.0  | 0.27 | 0.19 |
| 2.0 | 9.0  | 0.29 | 0.18 |
| ... | ...  | ...  | ...  |
| ... | ...  | ...  | ...  |


This file should be given as a binary format.

An example of LUT is `example/lut_860_2130.bin`.

You can see the details by a command:

    $ od -f example/lut_860_2130.bin | less

LUT should be calculated by radiative transfer models (ex. RSTAR).

Utils
---

* lut\_plot.rb

make a Nakajima-King-like plot of the look-up table.

It makes a plot like following:

![example](example/example.png)


References
---

* [Teruyuki Nakajima and Michael D. King, 1990: Determination of the Optical Thickness and Effective Particle Radius of Clouds from Reflected Solar Radiation Measurements. Part I: Theory. J. Atmos. Sci., 47, 1878–1893.](http://journals.ametsoc.org/doi/abs/10.1175/1520-0469(1990)047%3C1878%3ADOTOTA%3E2.0.CO%3B2)
* [Rodgers, C. D., 2000: Inverse Methods for Atmospheric Sounding-Theory and Practice.](http://www.worldscientific.com/worldscibooks/10.1142/3171)
* [Hiroshi Akima, 1970: A new method of interpolation and smooth curve fitting based on local procedures. JACM, 17.4, 589–602.](http://dl.acm.org/citation.cfm?id=321609)

