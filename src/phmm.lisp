;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;
;;;; PHMM definition - Pair Hidden Markov Models
;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :net.ashrentum.cl-hmm)

(eval-when (:compile-toplevel :load-toplevel :execute)

  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Model properties
  ;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  (deftype B-2streams-array ()
    `(prob-array (* * *)))

  );;;;end eval-when

;;;
;;; TODO create hmm-base class and make hmm-simple and phmm to inherit from it
;;;
(def-hmm-type phmm (hmm) *phmm-types* phmm-slots
  ((name
    :initarg :name :initform nil :accessor hmm-name)
   ;;;
   (S ;individual states
    :initarg :S :type vector-states :initform nil :accessor hmm-states)
   (S-hash ;state symbols -> index
    :initarg :S-hash :type hash-table :accessor hmm-states-hash)
   (N ;number of states
    :initarg :N :type cbook-state :accessor hmm-no-states)
   ;;;
   (L ;left stream alphabet
    :initarg :V :type alphabet :initform (error "Must set the left stream alphabet") :accessor hmm-alphabet-left)
   (L-hash ;left stream alphabet symbols -> index
    :initarg :V-hash :type hash-table :accessor hmm-alphabet-left-hash)
   (L-size ;left stream alphabet size
    :initarg :M :type cbook-symbol :accessor hmm-alphabet-left-size)
   ;;;
   (R ;right stream alphabet
    :initarg :V :type alphabet :initform (error "Must set the right stream alphabet") :accessor hmm-alphabet-right)
   (R-hash ;right stream alphabet symbols -> index
    :initarg :V-hash :type hash-table :accessor hmm-alphabet-right-hash)
   (R-size ;right stream alphabet size
    :initarg :M :type cbook-symbol :accessor hmm-alphabet-right-size)
   ;;;
   (groups ;group names
    :type vector-state-groups :accessor hmm-groups)
   (no-groups ;groups are the states which emission probs are tied
    :type cbook-state :initform 0 :accessor hmm-no-groups)
   (state-groups
    :type cbook-states :accessor hmm-state-groups)
   (PE ;initial state distribution, called PI in Rabiner //PE is PI in Phoenician
    :initarg :PE :type PE-vec :initform (error "Must set the initial probabilities") :accessor hmm-init)
   (A ;state transition probability distribution
    :initarg :A :type A-array :initform (error "Must set the transition probabilities") :accessor hmm-trans)
   (iA-from ;list, transitions from states to states
    :type itrans :accessor hmm-itrans-from)
   (iA-to ;list, transitions to states from states
    :type itrans :accessor hmm-itrans-to)
   (B ;left&right pair observation probability distribution
    :initarg :B :type B-2streams-array :initform (error "Must set the pair emission probabilities") :accessor hmm-emis)))


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Initialization
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod initialize-instance :after ((hmm phmm) &key)
  (phmm-slots (S S-hash N L L-hash L-size R R-hash R-size A) hmm
    ;;states hash, TODO define this in common place
    (do ((i 1 (1+ i)))
        ((> i N))
      (setf (gethash (state-name (aref S i)) S-hash) i))

    ;;L-hash, TODO define this in common place
    (do ((i 1 (1+ i)))
        ((> i L-size))
      (setf (gethash (aref L i) L-hash) i))

    ;;R-hash, TODO define this in common place
    (do ((i 1 (1+ i)))
        ((> i L-size))
      (setf (gethash (aref R i) R-hash) i))

    (multiple-value-bind (no-groups groups state-groups) (define-groups S)
                         (setf (hmm-no-groups hmm) no-groups)
                         (setf (hmm-groups hmm) groups)
                         (setf (hmm-state-groups hmm) state-groups))

    (setf (hmm-itrans-from hmm) (trans-array->itrans A))
    (setf (hmm-itrans-to hmm) (trans-array->itrans A t))
    (hmm-state-properties-set hmm)))
