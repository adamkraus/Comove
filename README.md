# Comove
A python package for finding comoving neighbors to a target star and returning useful information regarding youth (via FindFriends), as well as for quantitatively assessing the probability that a close-in neighbor is consistent with being a binary companion or a chance-aligned field star (via BinProb).<br/>

This package offers a quick reconnaissance of whether a potential young star has other young, comoving friends nearby and for assessing whether those friends are binary companions.<br/>

Written by Adam L. Kraus (UT Austin) and packaged by Aaron Rizzuto (UT Austin)<br/>

# What Comove/FindFriends does:
1. The code will query Gaia DR3 for stars within a volume radius equivalent to the defined search sphere, and with tangential velocities (from proper motion) within the input velocity difference range <br/>
2. Query 2MASS, GALEX, ROSAT, and WISE on the resulting neighbours <br/>
3. Make a series of useful plots displaying the search results <br/>
4. Output a table for all the neighbors containing the combined results <br/>

Currently only works with python 3<br/>

# DEPENDENCIES:
All available via pip:<br/>
galpy <br/>
astroquery<br/>
astropy<br/>
scipy<br/>
cartopy (and all of its dependencies)<br/>
matplotlib<br/>

# INSTALLATION:
To install the python 3 version:<br/>
git clone https://github.com/adamkraus/Comove.git<br/>
cd Comove<br/>
python setup.py build<br/>
python setup.py install<br/>

You can then run the code in the manner described in example.py<br/>

# NOTE

This package also contains the beta version of another function, BinProb, that considers the relative astrometry and photometry of a candidate binary companion and probabilistically assesses whether it is a bound binary companion or a chance-aligned field star. This code is under development, and it will be described in a future publication by A. Kraus, but users are welcome to experiment with it in the meatime. Caveat emptor, and remember, the warnings come after the spells.

# CITATION

If you use FindFriends, please cite Tofflemire et al. (2021), "TESS Hunt for Young and Maturing Exoplanets (THYME). V. A Sub-Neptune Transiting a Young Star in a Newly Discovered 250 Myr Association". This is where the motivation, functionality, and prototypical use-case are described. If you use BinProb, cite Deacon & Kraus (2020), as this is where much of the core mathematical formalism was first described.
