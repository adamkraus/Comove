# Comove
A python package for finding comoving neighbours to a target star and returning useful information regarding youth.<br/>
Great for seeing if a young looking star has other young, comoving friends nearby.<br/>

Written by Prof. Adam L. Kraus (UT Austin) and packaged by Aaron Rizzuto (UT Austin)<br/>

# What Comove does:
1. The code will query Gaia dr2 for stars within radius equivalent to the defined on-sky window, and within the input velocity difference range <br/>
2. Query 2MASS, Galex and Wise on the resulting neighbours <br/>
3. Make a series of useful plots displaying the search results <br/>
4. Output a table for all the neighbours containing the combined results <br/>

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
