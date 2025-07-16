import numpy as np
from astropy.io import fits
import os
import sys
########################################################################
def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")
########################################################################


reltrans_table = os.environ['RELTRANS_TABLES']
xillver_table_name = ['xillverCp_v3.4.fits', 'xillver-a-Ec5.fits', 'xillverD-5.fits']

print()
print('This script is going to create new version of the xillver tables')
print('in the same folder xillver tables (keep in mind the disk space)')
print()
query_yes_no('Do you agree?')

for table_name in xillver_table_name:    
    fits_image_filename = str(reltrans_table) + '/' + str(table_name)
    try:
        with fits.open(fits_image_filename) as hdul:
            spectra_extension = hdul['SPECTRA'].data
            length = len(spectra_extension.field(0))

            if (table_name == 'xillverCp_v3.4.fits'):
                logxi = spectra_extension.field(0)[:,2].reshape(length,1) # logxi ionisation values of the spectra
                logne = spectra_extension.field(0)[:,4].reshape(length,1) #  logne density values of the spectra
                spectra_extension['INTPSPEC'] = spectra_extension['INTPSPEC'] / 10**(logxi + logne - 15)
            if (table_name == 'xillver-a-Ec5.fits'):
                logxi = spectra_extension.field(0)[:,2].reshape(length,1) # logxi ionisation values of the spectra
                spectra_extension['INTPSPEC'] = spectra_extension['INTPSPEC'] / 10**(logxi)
            if (table_name == 'xillverD-5.fits'):
                logxi = spectra_extension.field(0)[:,2].reshape(length,1) # logxi ionisation values of the spectra
                logne = spectra_extension.field(0)[:,3].reshape(length,1) #  logne density values of the spectra
                spectra_extension['INTPSPEC'] = spectra_extension['INTPSPEC'] / 10**(logxi + logne - 15)

            norm_table_name = fits_image_filename[:-5] + '_normalised2.fits'
            hdul.writeto(norm_table_name)
        print(f'The script created a new table called {norm_table_name}')
        print()
    except:
        print(f'There is no table called {fits_image_filename}')

