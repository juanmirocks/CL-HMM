Copyright (C) Juan Miguel Cejuela <juanmi@jmcejuela.com>

--------------------------------------------------------------------------------
CL-HMM: HMM Library for Common Lisp
--------------------------------------------------------------------------------

* Compatibility: ANSI Common Lisp. Tested on SBCL 1.1.11 and Allegro 8.0 Free Express Edition
* Dependencies: [jmc.cl.utils system](https://github.com/jmcejuela/jmc.cl.utils)
* Development Times:
  * 2008 July-September
  * 2013-


Features:
--------------------------------------------------------------------------------

* Discrete observation densities.
* Alphabet symbols of any kind.
* Exponential state duration densities.
* Homogeneous HMMs.
* First order chains.
* Comparable efficiency to GHMM written in C (1x - 2x slower)
* HMMs:
  * Tied emission parameters
  * Finite and infinite HMMs
  * Forward/Backward scaled, Viterbi in log
  * Baum-Welch (scaled) for multiple labeled sequences and normalized noise

* PHMMs (Pair HMMs) -- **ongoing development**
  * Forward & backward pure algorithms (non-scaled)
  * Sequence translations: X --> Y
  * (soon) Baum-Welch pure
  * (soon) Viterbi
  * (soon) Scaled versions of the algorithms
  * (not soon; help anyone?) Tied emission parameters


Comments:
--------------------------------------------------------------------------------

- Followed Rabiner's notation in code. Otherwise properly indicated.


Installation:
--------------------------------------------------------------------------------

1. Download & load latest version of [jmc.cl.utils system](https://github.com/jmcejuela/jmc.cl.utils)
2. Download & load latest version of CL-HMM's code
3. Make ASDF locate the downloaded packaged, for instance: `(setf asdf:*central-registry* (append '(#p"...local path..." :*central-registry*)))`
4. `(asdf:operate 'asdf:load-op :cl-hmm)`
