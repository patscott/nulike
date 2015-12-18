nulike
======

Neutrino telescope likelihood tools


About
--

Author: Pat Scott (p.scott@imperial.ac.uk)  
Additional numerical routines contributed by Chris Savage (chris@savage.name)

Nulike is software for including full event-level information in
likelihood calculations for neutrino telescope searches for dark matter
annihilation.

Full details can be found in the papers:  
  1. Scott, Savage, Edsjö & IceCube Collaboration 2012, JCAP 11:057, [arXiv:1207.0810](http://arxiv.org/abs/1207.0810)
  2. IceCube Collaboration 2016, JCAP xx:xxx, [arXiv:1512.xxxxx](http://arxiv.org/abs/1512.xxxxx)
  
If you use nulike for preparing a paper, please cite both of these.


Compilation
--

To compile nulike as a static library, edit [makefile] and do

  make

To compile nulike as a shared library, instead do

  make libnulike.so

These libraries can then be linked against whatever driver program you
want to use to provide the neutrino flux at Earth from your favourite
dark matter model.  I recommend DarkSUSY 5.1.3 (www.darksusy.org) or
later, but in principle anything that computes a neutrino flux is OK
(just make sure it is threadsafe if you want to run nulike
multi-threaded).


Testing
--

To compile the nulike test programs, you must first install DarkSUSY and
enter the path to its library in [makefile].  Then, depending on which
demo you want to run, you can do

  make nulike_test  
  make nulike_test_wimp  
  make nulike_test_mssm25

Please refer to the source code of these examples for descriptions of
options, etc.  If the output from nulike_test_wimp seems to match the
contents of [nulike_test_wimp.diffme], then things are working properly.


Generating partial likelihoods
--

To (re-)compute partial likelihoods for later reuse in fast likelihood
calculations, you must first install nusigma (found from the WIMPSim
page at http://copsosx03.fysik.su.se/wimpsim/code.html) and enter the
path to its library in [makefile].  Then,

  make nulike_prep

and

  ./nulike_prep

to get full usage instructions.
