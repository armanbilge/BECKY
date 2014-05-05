# BECKY
## Bayesian Estimation of Coevolutionary KrYteria

A plugin for [BEAST1](//beast-mcmc.googlecode.com/) that enables inference of cophylogenies and, in particular, coevolutionary phenomenon in an ecological context.

[**DOWNLOAD**](//github.com/coevolution/BECKY/releases)

[![Build Status](//travis-ci.org/armanbilge/BECKY.svg?branch=master)](//travis-ci.org/armanbilge/BECKY)

## Installation and Usage

To install BECKY, simply create a `plugins/` folder in the same location as your BEAST XML file and place `org.ithinktree.becky.BECKY.jar` in it.

1. To setup an analysis, first use BEAUti to setup two BEAST XMLs, one for the host organism and another for the symbiont.

2. Next, prepare an associations file in the following format:
   ```
   symbiont1.1 <tab> host1
   symbiont1.2 <tab> host1
   symbiont2.1 <tab> host2
   ```

3. Run the included python script `SetupCophylogenyAnalysis.py` with the following arguments:
   ```
   python SetupCophylogenyAnalysis.py -a <associations-file> -j <host-prefix> -k <host-xml>
                                       -s <symbiont-prefix> -t <symbiont-xml> > <output-xml>
   ```
   **Note that you must redirect `stdout` to a new XML file.**
   
   For example:
   ```
   python SetupCophylogenyAnalysis.py -a assoc.txt -j gopher -k gopher.xml
                                       -s louse -t louse.xml > coevolution.xml
   ```
   

4. BEAST can be started normally. For BEAST to succesfully find BECKY, make sure to place the jar file in a `plugins/` folder at the same location as your input XML.
