.. Environ documentation 'on handling output files', created
   by Matthew Truscott on Mon Apr 8 2019.

Output Files: Reading and Automation
====================================

Unlike Quantum-ESPRESSO, which parses the entire output into an XML file, Environ currently relies solely on
the output file and an associated environ-debug file (if the appropriate verbosity is set in the input).
The reason for this is that Environ is currently a plugin and does not overwrite or modify any of the base
Quantum-ESPRESSO files. Hence any additional output is only currently printed out into the terminal (or piped
into an output file), or in environ specific output files. This may be dealt with more elegantly in the future,
but for now, when using Environ, it is more appropriate to rely on the direct output of the program. One can
feed this into a file and parse that file, typically using grep on certain keywords, or using scripting.
Users have developed post-processing tools for extracting information for personal use. These can be 
particularly attractive when running many simulations for a particular investigation and as such, it is likely
that we will release some helpful post-processing tools along with future releases. 


