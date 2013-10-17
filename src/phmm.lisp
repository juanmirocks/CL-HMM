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
   (N ;number of states
    :initarg :N :type cbook-state :accessor hmm-no-states)
   (S-hash ;state symbols -> index
    :initarg :S-hash :type hash-table :accessor hmm-states-hash)
   ;;;
   (L ;left stream alphabet
    :initarg :L :type alphabet :initform (error "Must set the left stream alphabet") :accessor hmm-alphabet-left)
   (L-size ;left stream alphabet size
    :initarg :L-size :type cbook-symbol :accessor hmm-alphabet-left-size)
   (L-hash ;left stream alphabet symbols -> index
    :initarg :L-hash :type hash-table :accessor hmm-alphabet-left-hash)
   ;;;
   (R ;right stream alphabet
    :initarg :R :type alphabet :initform (error "Must set the right stream alphabet") :accessor hmm-alphabet-right)
   (R-size ;right stream alphabet size
    :initarg :R-size :type cbook-symbol :accessor hmm-alphabet-right-size)
   (R-hash ;right stream alphabet symbols -> index
    :initarg :R-hash :type hash-table :accessor hmm-alphabet-right-hash)
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

(defmethod print-object ((object phmm) stream)
  (declare (stream stream))
  (labels ((print-A (A S)
             (let ((out ""))
               (dolist-index (trans (multiple-value-bind (a b) (trans-stats A S)
                                      (setf out (format nil " (~a): ~&" a)) b) i out)
                 (setf out
                       (concatenate 'string out (format nil "~5T~a (~a):  ~{~a:~3$~^  ~}~%"
                                                        (car trans) (second trans) (cddr trans)))))))
           (print-B ()
             "...")) ;;TODO

  (print-unreadable-object (object stream :type t)
    (phmm-slots (name S N L L-size R R-size PE A) object
      (let ((no-begins) (begins))
        (multiple-value-setq (no-begins begins) (init-stats PE S))
        (format stream "~a (~d,~d,~d)~%  S: ~a~2%  L: ~a~%  R: ~a~2%  PE (~a): ~{~a:~3$~^  ~}~2%  A~a~%  B: ~a~%"
                name N L-size R-size
                S
                L
                R
                no-begins
                begins
                (print-A A S)
                (print-B)))))))


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Initialization
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod initialize-instance :after ((hmm phmm) &key)
  (phmm-slots (S S-hash N L L-hash L-size R R-hash R-size A) hmm

    ;;states hash, TODO define this in common place
    (do* ((i 0 i+1)
          (i+1 1 (1+ i+1)))
         ((= i N))
      (setf (gethash (state-name (aref S i)) S-hash) i+1))

    ;; ;;L-hash, TODO define this in common place
    (do* ((i 0 i+1)
          (i+1 1 (1+ i+1)))
        ((= i L-size))
      (setf (gethash (aref L i) L-hash) i+1))

    ;; ;;R-hash, TODO define this in common place
    (do* ((i 0 i+1)
          (i+1 1 (1+ i+1)))
         ((= i R-size))
      (setf (gethash (aref R i) R-hash) i+1))

    (multiple-value-bind (no-groups groups state-groups) (define-groups S)
                         (setf (hmm-no-groups hmm) no-groups)
                         (setf (hmm-groups hmm) groups)
                         (setf (hmm-state-groups hmm) state-groups))

    (setf (hmm-itrans-from hmm) (trans-array->itrans A))
    (setf (hmm-itrans-to hmm) (trans-array->itrans A t))
    ;;(hmm-state-properties-set hmm) TODO
    ))

(defun make-phmm (N L-list R-list model &key name (L-alphabet-type T) (R-alphabet-type T) (model-spec :complete))
  "Make a phmm. Two ways to specify the model parameters as follows:
  N:
  L-size:
  R-size:
  alphabet: list of emissions symbols (eg, '(A C G T))
  model: depending on model-spec this, list has 2 different forms
  model form 1: (model-spec = :complete)
    states: state names and labels
    init: initial probs
    trans: trans probs
    emis: trans probs
    example: '(((:fair #\F) (:biased #\B)) (1 0) ((.95 .05) (.15 .85)) ((1/6 ...) (1/2 1/10 ...)))

  model form 2: (model-spec = :relevant)
    NOT IMPLEMENTED"

  (let* ((S (make-array N :element-type 'state :initial-element (list (make-typed-array 8 'bit 0))))
         (states)
         (S-hash)
         ;;
         (L-size (length L-list))
         (L (make-array L-size :element-type L-alphabet-type))
         (L-hash)
         ;;
         (R-size (length R-list))
         (R (make-array R-size :element-type R-alphabet-type))
         (R-hash)
         ;;
         (PE (make-typed-array N 'prob-float +0-prob+))
         (A (make-typed-array (list N N) 'prob-float +0-prob+))
         (B (make-typed-array (list N (1+ L-size) (1+ R-size)) 'prob-float +0-prob+)))
;;;depends on model specification
    (ecase model-spec
      (:complete
       (setf states (first model))
       (setf S-hash (make-hash-table-with-list states :test 'equal)) ;S-hash
       (dolist-index (a (second model) i) (setf (aref PE i) (prob a))) ;init probs
       (dolist-index (state-transitions (third model) i)
         (dolist-index (at state-transitions j)
           (setf (aref A i j) (prob at))))
       (dolist-index (state-emissions (fourth model) i)
         (dolist-index (left state-emissions j)
           (dolist-index (p left z)
             (setf (aref B i j z) (prob p))))))
      (:relevant
       (error "relevant specification: NOT IMPLEMENTED")))
;;;same for all model specifications
    (dolist-index (symbol L-list i) (setf (aref L i) symbol))
    (dolist-index (symbol R-list i) (setf (aref R i) symbol))
    (setf L-hash (make-hash-table-with-list L-list)) ;TODO watch out with the indexes, must start in 1
    (setf R-hash (make-hash-table-with-list R-list)) ;TODO watch out with the indexes, must start in 1

    (dolist-index (e states i) ;the states
      (setf (aref S i) (cons (make-typed-array 8 'bit 0) (if (listp e) e (list e +label-null+)))))
    (!normalize-vector PE) ;normalize probs
    (!normalize-2dmatrix-by-row A)
    ;;emission probability(epsilon & epsilon) == 0
    (dotimes (i N)
      (setf (aref B i 0 0) +0-prob+))
    (!normalize-3dmatrix-by-row B)

    (make-instance 'phmm :name name
                   :S S :N N :S-hash S-hash
                   :L L :L-size L-size :L-hash L-hash
                   :R R :R-size R-size :R-hash R-hash
                   :PE PE :A A :B B)))

(defun make-random-phmm
    (N L-size R-size &key (eccentricity 2) (states (range N)) (L-list (range L-size)) (R-list (range R-size)) name (L-alphabet-type T) (R-alphabet-type T))
  "Make a phmm with no biased info
  ...
  eccentricity: real, eccentricity for randomly-generated probabilities. The bigger the more dispair. If 0, the probs. are uniform
  states: list / optional list of states info (eg, ((:fair #\F) (:biased #\B)))
  ..."
  (make-phmm N L-list R-list
             (list
              states
              (make-list-meval N (expt (random 1.0) eccentricity))
              (make-list-meval N (make-list-meval N (expt (random 1.0) eccentricity)))
              (make-list-meval N (make-list-meval (1+ L-size) (make-list-meval (1+ R-size) (expt (random 1.0) eccentricity)))))
             :name name :L-alphabet-type L-alphabet-type :R-alphabet-type R-alphabet-type :model-spec :complete))


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Common Methods
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod hmm-copy ((hmm phmm))
  (phmm-slots (S N L L-size R R-size) hmm
    (make-instance (type-of hmm)
                   :S S
                   :N N
                   :S-hash (make-hash-table :size N :test 'equalp)

                   :L L
                   :L-size L-size
                   :L-hash (make-hash-table :size L-size :test 'equalp)

                   :R R
                   :R-size R-size
                   :R-hash (make-hash-table :size R-size :test 'equalp)

                   :PE (copy-seq (hmm-init hmm))
                   :A (copy-matrix (hmm-trans hmm))
                   :B (copy-array (hmm-emis hmm)))))


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Forward & Backward
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod forward ((hmm phmm) obs-c)
  "
@param hmm: pair hidden markov model
@param obs-c: cbook-encoded pair observation, list of 2 elements

@return (1) probability of observation pair given hmm model
@return (2) generated alpha 3d matrix
"
  (declare (optimize (speed 3) (safety 0)))
  (phmm-slots (N PE A B iA-to) hmm
    (let* ((x (first obs-c))
           (y (second obs-c))
           (n (length x))
           (m (length y))
           (alpha (make-typed-array `(,N (1+ ,n) (1+ ,m)) 'prob-float +0-prob+)))
      (declare (fixnum n m)
               (cbook-alphabet x y))

      ;; Initialisation
      ;; -------------------------------------------------------------------------
      (loop for j below N do
           (setf (aref alpha j 1 0) (prob (* (aref PE j) (aref B j (aref x 0) +epsilon-cbook-index+))))
           (setf (aref alpha j 0 1) (prob (* (aref PE j) (aref B j +epsilon-cbook-index+ (aref y 0)))))
           (setf (aref alpha j 1 1)
                 (prob (+
                        (* (aref PE j) (aref B j (aref x 0) (aref y 0)))
                        (* (loop for i below N sum (* (aref A i j) (aref alpha i 0 1))) (aref B j (aref x 0) +epsilon-cbook-index+))
                        (* (loop for i below N sum (* (aref A i j) (aref alpha i 1 0))) (aref B j +epsilon-cbook-index+ (aref y 0)))))))


      ;; Induction
      ;; -------------------------------------------------------------------------
      ;; TODO gain speed by checking whether 0 == b_j(x_l, y_r)
      (loop for j below N do
           (loop for l from 0 to n do
                (loop for r from 0 to m do
                     (when (<= 2 (min l r))
                       (setf (aref alpha j l r)
                             (prob
                              (+
                               (* (loop for i below N sum (* (aref A i j) (aref alpha i (1- l) r))) (aref B j (aref x l) (aref y r)))
                               (* (loop for i below N sum (* (aref A i j) (aref alpha i (1- l) r))) (aref B j (aref x l) +epsilon-cbook-index+))
                               (* (loop for i below N sum (* (aref A i j) (aref alpha i l (1- r)))) (aref B j +epsilon-cbook-index+ (aref y r))))))))))


      ;; Termination
      ;; -------------------------------------------------------------------------
      (values
       (loop for j below N sum (aref alpha n m))
       alpha))))
