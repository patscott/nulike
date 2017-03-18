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
  2. IceCube Collaboration 2016, JCAP 04:022, [arXiv:1601.00653](http://arxiv.org/abs/1601.00653)

If you use nulike for preparing a paper, please cite both of these.  If you use the 79-string IceCube likelihoods, don't forget to also cite the original IC79 WIMP paper:
  3. IceCube Collaboration 2013, PRL 110:131302, [arXiv:1212.4097](http://arxiv.org/abs/1212.4097)


Compilation
--

To compile nulike as a static library, edit [the makefile](makefile) and do
```
  make
```

To compile nulike as a shared library, instead do
```
  make libnulike.so
```

These libraries can then be linked against whatever driver program you
want to use to provide the neutrino flux at Earth from your favourite
dark matter model.  I recommend DarkSUSY 5.1.3 (www.darksusy.org) or
later, but in principle anything that computes a neutrino flux is OK
(just make sure it is threadsafe if you want to run nulike
multi-threaded).


Testing
--

To compile the nulike test programs, you must first install DarkSUSY and
enter the path to its library in [the makefile](makefile).  Then, depending on which
demo you want to run, you can do
```
  make nulike_test
  make nulike_test_wimp
  make nulike_test_mssm25
```
Please refer to the source code of these examples for descriptions of
options, etc.  If the output from nulike_test_wimp seems to match the
contents of [tests/nulike_test_wimp.diffme](tests/nulike_test_wimp.diffme), then things are working
properly.


Generating partial likelihoods
--

To (re-)compute partial likelihoods for later reuse in fast likelihood
calculations, you must first install nusigma (found from the WIMPSim
page at http://copsosx03.fysik.su.se/wimpsim/code.html) and enter the
path to its library in [the makefile](makefile).  Then,
```
  make nulike_prep
```
and
```
  ./nulike_prep
```
to get full usage instructions.



License
--
Copyright (c) 2015, Patrick Scott and Christopher Savage
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
