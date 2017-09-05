# Author: Willis O'Leary

def warn(desc):
    """Prints a warning."""
    print 'WARNING: {0}'.format(desc)
def error(desc):
    """Prints an error and exits."""
    print 'ERROR: {0}'.format(desc)
    exit()

""" Atom masses."""
masses = {'H':1.008,  'He':4.003, 'Li':6.941,  'Be':9.012182,
          'B':10.811 , 'C':12.0107, 'N':14.00674,'O':15.9994, 
          'F':18.9984032, 'Ne':20.1797,'Na':22.98977,'Mg':24.3050,
          'Al':26.981538, 'Si':28.0855,'P':30.973761,'S':32.066, 
          'Cl':35.4527, 'Ar':39.948,'K':39.0983,'Ca':40.078, 
          'Pd':106.42, 'Mo': 95.94, 'Nb':92.91, 'Te':127.6, 
          'V':50.94, 'Ir':192.2, 'Cu':63.546}

kb = 8.3148e-7 # Boltzmann constant in  angst^2 amu / fs^2 K






