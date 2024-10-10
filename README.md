# orca_cp_writer
Generates an Orca counterpoise interaction energy compound job input

This has been tested with Orca version 5, and applies the Boys-Bernardi counterpoise correction for a single point energy as detailed in [Mol. Phys. 19, 553-566 (1970)](https://doi.org/10.1080/00268977000101561).

Assuming the file has been made executable, basic usage information is available via `orca_cp_write.py -h`

It will attempt to detect if an empirical dispersion correction has been requested and add the extra input lines required.
