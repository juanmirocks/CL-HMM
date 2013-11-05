Copyright (C) Juan Miguel Cejuela <juanmi@jmcejuela.com>

--------------------------------------------------------------------------------
CL-HMM: HMM Library for Common Lisp
--------------------------------------------------------------------------------

* Compatibility: ANSI Common Lisp. Tested on SBCL 1.1.11 and Allegro 8.0 Free Express Edition
* Dependencies: [jmc.cl.utils system](https://github.com/jmcejuela/jmc.cl.utils)
* Development Times: 2008 July-September, 2013-


Features:
--------------------------------------------------------------------------------

* Discrete observation densities.
* Alphabet symbols of any kind.
* Exponential state duration densities.
* Homogeneous HMMs.
* First order chains.

* HMMs:
  * Tied emission parameters
  * Finite and infinite HMMs
  * Forward/Backward scaled, Viterbi in log
  * Baum-Welch (scaled) for multiple labeled sequences and normalized noise
  * Comparable efficiency to GHMM written in C (1x - 2x slower)

* PHMMs (Pair HMMs) -- **ongoing development**
  * Forward & backward pure algorithms (non-scaled)
  * (ONGOING) Baum-Welch pure
  * (ONGOING) Sequence translations: X -> Y


Comments:
--------------------------------------------------------------------------------

- Followed Rabiner's notation in code. Otherwise properly indicated.
