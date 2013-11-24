Copyright (C) Juan Miguel Cejuela <juanmi@jmcejuela.com>

--------------------------------------------------------------------------------
CL-HMM: HMM Library for Common Lisp
--------------------------------------------------------------------------------

* Compatibility: ANSI Common Lisp. Tested on SBCL 1.1.11
* Dependencies: [jmc.cl.utils system](https://github.com/jmcejuela/jmc.cl.utils)
* Development Times:
  * 2008 July-September
  * 2013-


Features:
--------------------------------------------------------------------------------

* Discrete observation densities
* Alphabet symbols of any kind
* Exponential state duration densities
* Homogeneous HMMs
* First order chains
* Comparable efficiency to GHMM written in C (1x - 2x slower)
* HMMs:
  * Forward & backward in probability space and scaled
  * Viterbi in probability and log space
  * Baum-Welch training (probability and scaled) for multiple sequences (optionally labeled)
  * Tied emission parameters
  * Finite and infinite HMMs

* PHMMs (Pair HMMs): **ongoing development**
  * Forward & backward in probability and log space
  * (soon) Viterbi in probability and log space
  * (soon) Baum-Welch training (probability and log) for multiple sequences
  * Sequence translations: given input X --> generate output Y


Comments:
--------------------------------------------------------------------------------

- Rabiner's notation is followed in code. Otherwise properly indicated.


Installation:
--------------------------------------------------------------------------------

1. You must have an installed Common Lisp implementation. For instance: [SBCL](http://www.sbcl.org/)
2. Download latest version of [jmc.cl.utils system](https://github.com/jmcejuela/jmc.cl.utils)
3. Download latest version of CL-HMM's code
4. Make ASDF locate the downloaded packaged, for instance: `(setf asdf:*central-registry* (append '(#p"...local path..." :*central-registry*)))`
5. Load & compile the package: `(asdf:operate 'asdf:load-op :cl-hmm)`
