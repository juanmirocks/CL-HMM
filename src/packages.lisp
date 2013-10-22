;; Author: Juan Miguel Cejuela
;; Created: Wed Jul  9 19:13:05 2008 (CEST)
;; Version:
;; Last-Updated: 2011-08-11
;;           By: Juan Miguel Cejuela
;;     Update #: 9

(defpackage :net.ashrentum.cl-hmm
  (:use :cl :jmcejuela)
  (:nicknames :cl-hmm)
  (:export

   ;;;cl-hmm
   :hmm
   :hmm-copy
   :hmm-correctp
   :hmm-compatiblep
   :hmm-run

   ;;;hmm-simple
   :hmm-simple
   :hmm-infinite
   :hmm-finite

   :hmm-state-labels

   :hmm-name
   :hmm-no-states
   :hmm-no-emissions
   :hmm-alphabet
   :hmm-alphabet-hash
   :hmm-states
   :hmm-states-hash
   :hmm-no-groups
   :hmm-groups
   :hmm-state-groups
   :hmm-init
   :hmm-trans
   :hmm-itrans-from
   :hmm-itrans-to
   :hmm-emis

   :make-hmm-simple
   :make-random-hmm-simple
   :get-statename-index
   :begin-state-3p
   :begin-state-2p
   :begin-state-p
   :end-state-3p
   :end-state-2p
   :end-state-p
   :silent-state-3p
   :silent-state-2p
   :silent-state-p
   :invalid-state-3p
   :invalid-state-2p
   :invalid-state-p
   :tied-state-3p
   :tied-state-2p
   :tied-state-p
   :fixed-state-3p
   :fixed-state-2p
   :fixed-state-p

   :cbook-encode
   :cbook-decode

   ;;;phmm
   :phmm

   :make-phmm

   ;;forward&backward
   :forward
   :likelihood
   :forward-scl
   :loglikelihood
   :backward
   :backward-scl
   :posterior-probs-labels
   ;;viterbi
   :viterbi
   :viterbi-log

   ;;baum-welch
   :baum-welch
   :baum-welch-scl
   :hmm-estimate-bws

   ;;training
   :soften-label-boundaries

   ;;alphabets
   :*protein-code*
   :*adn-code*
   :*arn-code*

   ;;hmm-files
   :hmm-save
   :hmm-open))
