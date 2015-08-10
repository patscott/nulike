#!/bin/python
# Makes SLHA files into MSSM-25 models in DarkSUSY format.
# Not at all SLHA-compliant, smart or fast; just quick to write.
# Pat Scott

import re
import sys  
import collections

# No empties from re.split
def neatsplit(regex,string):
    return [x for x in re.split(regex,string) if x != '']

# Set up the output order as 
#.....am1...... .....am2...... .....am3......
#......amu..... .....ama...... ....atanbe.... 
#....amsqL1.... ....amsql2.... ....amsql3.... ....amsqRu.... ....amsqRc.... ....amsqRt.... ....amsqRd.... ....amsqRs.... ....amsqRb....
#....amslL1.... ....amslL2.... ....amslL3.... ....amslRe.... ...amslRmu.... ...amslRtau...
#.....atm...... .....abm...... ....ataum..... ....aemum.....
entries = collections.OrderedDict() 
entries['M_1(Q)'] = 0.0 
entries['M_2(Q)'] = 0.0
entries['M_3(Q)'] = 0.0
entries['mu(Q)MSSM'] = 0.0
entries['mA^2(Q)'] = 0.0
entries['tan'] = 0.0
entries['M_q1L'] = 0.0
entries['M_q2L'] = 0.0
entries['M_q3L'] = 0.0
entries['M_uR'] = 0.0
entries['M_cR'] = 0.0
entries['M_tR'] = 0.0
entries['M_dR'] = 0.0
entries['M_sR'] = 0.0
entries['M_bR'] = 0.0 
entries['M_eL'] = 0.0
entries['M_muL'] = 0.0
entries['M_tauL'] = 0.0
entries['M_eR'] = 0.0
entries['M_muR'] = 0.0
entries['M_tauR'] = 0.0
entries['A_t(Q)'] = 0.0
entries['A_b(Q)'] = 0.0
entries['A_tau(Q)'] = 0.0
entries['A_mu(Q)'] = 0.0

def main(argv):
  f = open(argv[0], 'r')
  for line in f:
    splitline = neatsplit('\s|\#',line)
    if len(splitline) < 3 or splitline[0] == 'BLOCK': continue
    for key in entries:
      if key in splitline[2]:
        if key == 'mA^2(Q)':
          entries[key] = pow(float(splitline[1]),0.5)
        else:
          entries[key] = float(splitline[1])
      if len(splitline) < 4: continue
      if key in splitline[3]: 
        entries[key] = float(splitline[2])
      elif key in splitline[2]+' '+splitline[3]:
        entries[key] = float(splitline[1])
  f.close()
  print "{:12s}".format(" "+argv[1]),
  for val in entries.values():
   print "{:14.5e}".format(val),  
   
# Handle command line arguments
if __name__ == "__main__":
   main(sys.argv[1:])
