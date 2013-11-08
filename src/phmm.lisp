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

(defun make-uniform-phmm
    (N L-size R-size &key (states (range N)) (L-list (range L-size)) (R-list (range R-size)) name (L-alphabet-type T) (R-alphabet-type T))
  "Same as make-random-hmm-simple but with an uniform distribution for all probabilities, eccentricity 0"
  (make-random-phmm N L-size R-size :eccentricity 0 :states states :L-list L-list :R-list R-list :name name :L-alphabet-type L-alphabet-type :R-alphabet-type R-alphabet-type))


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Common Methods
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod hmm-correctp ((hmm phmm))
  (warn "TODO, Not fully implemented")
  (let ((buffer (make-array +stream-correctness-size+ :element-type 'character :adjustable t :fill-pointer 0))
        (correct t))

    (labels ((existsNegative (matrix)
               (dotimes (i (array-total-size matrix) nil)
                 (if (< (row-major-aref matrix i) 0) (return t))))
             (assume (condition errorMsg)
               (unless condition (setq correct nil) (format buffer errorMsg))))


      (phmm-slots (N PE A B) hmm
        (assume (not (existsNegative PE)) "Negative value in PE")
        (assume (not (existsNegative A)) "Negative value in A")
        (assume (not (existsNegative B)) "Negative value in B")
        (loop for i below N do
             (assume (zerop (aref B i 0 0)) "Non-zero! b_~i(epsilon, epsilon)"))))

    (values correct buffer)))

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

;;;this definition can change
(defmethod hmm-complexity ((hmm phmm))
  (declare (inline hmm-complexity))
  (phmm-slots (N L-size R-size) hmm
    (* N L-size R-size)))

(defmethod cbook-encode-left ((hmm phmm) x)
  "cbook-encode left string, that is, the input of a pair observation"
  (declare (optimize (speed 3)) ((vector) x))
  (let* ((size_x (the fixnum (length x)))
         (out_x (make-array size_x :element-type 'cbook-symbol)))
    (phmm-slots (L-hash) hmm
      (dotimes (i size_x out_x)
        (setf (aref out_x i) (gethash (aref x i) L-hash))))))

(defmethod cbook-encode-right ((hmm phmm) y)
  "cbook-encode right string, that is, the output of a pair observation"
  (declare (optimize (speed 3)) ((vector) y))
  (let* ((size_y (the fixnum (length y)))
         (out_y (make-array size_y :element-type 'cbook-symbol)))
    (phmm-slots (R-hash) hmm
      (dotimes (i size_y out_y)
        (setf (aref out_y i) (gethash (aref y i) R-hash))))))

(defmethod cbook-encode ((hmm phmm) observation)
  "cbook-encode pair observation"
  (declare (optimize (speed 3)))
  (let* ((x (first observation))
         (y (second observation)))
    (list (cbook-encode-left hmm x) (cbook-encode-right hmm y))))

(defmethod cbook-decode ((hmm phmm) observation)
  "cbook-encode pair observation"
  (declare (optimize (speed 3) (safety 0)))
  (let* ((in_x (first observation))
         (in_y (second observation))
         (size_x (the fixnum (length in_x)))
         (size_y (the fixnum (length in_y))))
    (declare (cbook-alphabet in_x in_y))
    (phmm-slots (L R) hmm
      (let ((out_x (make-array size_x :element-type (array-element-type L)))
            (out_y (make-array size_y :element-type (array-element-type R))))

        (list
         (dotimes (i size_x out_x)
           (setf (aref out_x i) (aref L (1- (aref in_x i))))) ;-1 since in L&R epsilon is not accounted for
         (dotimes (i size_y out_y)
           (setf (aref out_y i) (aref R (1- (aref in_y i))))))))))

(defun cbelt1 (seq i)
  "1-indexed cbook-encoded input sequence. If 0, return epsilon's index"
  (declare (inline cbelt1))
  (if (= i 0)
      +epsilon-cbook-index+
      (elt seq (1- i))))

(defmacro arefalpha (matrix dim1 dim2 dim3)
  "Alpha matrix accessor"
  `(if (or (< ,dim2 0) (< ,dim3 0))
       +0-prob+
       (aref ,matrix ,dim1 ,dim2 ,dim3)))

(defmethod !hmm-noisify ((hmm phmm) noise)
  (unless (zerop noise)
    (let ((confidence (coerce (- 1 noise) 'prob-float)))
      (phmm-slots (PE A B) hmm
        (!normalize-vector
         (!combine-float-arrays PE (!normalize-vector (make-random-array (array-dimensions PE) +1-prob+)) confidence t))
        (!normalize-2dmatrix-by-row
         (!combine-float-arrays A (!normalize-2dmatrix-by-row (make-random-array (array-dimensions A) +1-prob+)) confidence t))
        (!normalize-3dmatrix-by-row
         (!combine-float-arrays B (!normalize-3dmatrix-by-row (make-random-array (array-dimensions B) +1-prob+)) confidence t)))))
  hmm)

(defmethod hmm-save ((hmm phmm) filename &optional model-spec)
  (when (not (eql model-spec :complete)) (error "Only the model-spec ':complete' is supported"))
  (phmm-slots (name N L R S PE A B) hmm
    (let ((states (loop for i below N collect (list (state-name (aref S i)) (state-label (aref S i))))))
      (with-open-file (stream filename :direction :output :if-does-not-exist :create :if-exists :rename)
        (prin1
         `(make-phmm
           ,N
           ',(array->list L)
           ',(array->list R)
           ,`'(,states
               ,(array->list PE)
               ,(array->list A)
               ,(array->list B))
           :name ,name
           :L-alphabet-type ',(array-element-type L)
           :R-alphabet-type ',(array-element-type R)
           :model-spec ,model-spec)
         stream)))))

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
           (size_x (length x))
           (size_y (length y))
           (alpha (make-typed-array `(,N ,(1+ size_x) ,(1+ size_y)) 'prob-float +0-prob+)))
      (declare (fixnum size_x size_y)
               (cbook-alphabet x y))

      ;; Initialization
      ;; -------------------------------------------------------------------------
      (loop for j below N do
           (when (> size_x 0) (setf (aref alpha j 1 0) (prob (* (aref PE j) (aref B j (cbelt1 x 1) +epsilon-cbook-index+)))))
           (when (> size_y 0) (setf (aref alpha j 0 1) (prob (* (aref PE j) (aref B j +epsilon-cbook-index+ (cbelt1 y 1))))))
           (when (and (> size_x 0) (> size_y 0))
             (setf (aref alpha j 1 1)
                   (prob (+
                          (* (aref PE j) (aref B j (aref x 0) (aref y 0)))
                          (* (loop for i in (aref iA-to j) sum (* (aref A i j) (aref alpha i 0 1))) (aref B j (cbelt1 x 1) +epsilon-cbook-index+))
                          (* (loop for i in (aref iA-to j) sum (* (aref A i j) (aref alpha i 1 0))) (aref B j +epsilon-cbook-index+ (cbelt1 y 1))))))))


      ;; Induction
      ;; -------------------------------------------------------------------------
      ;; TODO gain speed by checking whether 0 == b_j(x_l, y_r)
      (loop for l from 0 to size_x do
           (loop for r from 0 to size_y do
                (when (<= 2 (max l r))
                  (loop for j below N do
                       (setf (aref alpha j l r)
                             (prob
                              (+
                               (* (loop for i in (aref iA-to j) sum (* (aref A i j) (arefalpha alpha i (1- l) (1- r)))) (aref B j (cbelt1 x l) (cbelt1 y r)))
                               (* (loop for i in (aref iA-to j) sum (* (aref A i j) (arefalpha alpha i (1- l) r))) (aref B j (cbelt1 x l) +epsilon-cbook-index+))
                               (* (loop for i in (aref iA-to j) sum (* (aref A i j) (arefalpha alpha i l (1- r)))) (aref B j +epsilon-cbook-index+ (cbelt1 y r))))))))))


      ;; Termination
      ;; -------------------------------------------------------------------------
      (values
       (loop for j below N sum (aref alpha j size_x size_y))
       alpha))))


(defmethod backward ((hmm phmm) obs-c)
  "
@param hmm: pair hidden markov model
@param obs-c: cbook-encoded pair observation, list of 2 elements

@return generated beta 3d matrix
"
  (declare (optimize (speed 3) (safety 0)))
  (phmm-slots (N A B iA-from) hmm
    (let* ((x (first obs-c))
           (y (second obs-c))
           (size_x (length x))
           (size_y (length y))
           (beta (make-typed-array `(,N ,(1+ size_x) ,(1+ size_y)) 'prob-float +0-prob+)))
      (declare (fixnum size_x size_y)
               (cbook-alphabet x y))

      (labels ((arefmat (matrix dim1 dim2 dim3)
                 (if (or (> dim2 size_x) (> dim3 size_y))
                     0
                     (aref matrix dim1 dim2 dim3)))
               (elt1 (seq i)
                 "1-indexed cbook-encoded input sequence. if i >= length(seq), return epsilon's index"
                 (if (> i (length seq))
                     +epsilon-cbook-index+
                     (elt seq (1- i)))))

        ;;Initialization
        ;; -------------------------------------------------------------------------
        (loop for i below N do
             (setf (aref beta i size_x size_y) (prob 1)))

        ;;Induction
        ;; -------------------------------------------------------------------------
        (loop for l from size_x downto 0 do
             (loop for r from size_y downto 0 do
                  (when (and (<= 1 (max l r)) (not (and (= l size_x) (= r size_y))))
                    (loop for i below N do
                         (setf (aref beta i l r)
                               (prob
                                (loop for j in (aref iA-from i) sum
                                     (* (aref A i j)
                                        (+ (* (arefmat beta j (1+ l) (1+ r)) (aref B j (elt1 x (1+ l))       (elt1 y (1+ r))))
                                           (* (arefmat beta j (1+ l) r     ) (aref B j (elt1 x (1+ l))       +epsilon-cbook-index+))
                                           (* (arefmat beta j l      (1+ r)) (aref B j +epsilon-cbook-index+ (elt1 y (1+ r))))))))))))))

      beta)))



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Translations
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod hmm-translate ((phmm phmm) X)
  "Translate X input sequence (left stream) -> Y output sequence (right stream)

  The Y translation is non-deterministic and is constructed following a random path
  determined by the transition and emission distributions of the model. Some randomly-generated
  paths may lead to dead ends, that is, where the X is not completely generated.

  @param phmm: model with probability distributions to use for translation
  @param X: cbook-encoded left-stream input sequence

@return(0) Y: a translation of X, cbook-encoded right-stream input sequence"

  (declare (cbook-alphabet X))

  (phmm-slots (PE A B) phmm
    (let ((size_x (length X))
          (init_state (select-random (accum-array PE 1 prob-float)))
          (accumA (accum-array A 2 prob-float))
          (accumB (accum-array B 3 prob-float)))
      (labels ((zerop-emission-prob (state l) (not (select-random accumB :indices-1 (list state l))))
               (rec (l state Y)
                 (if (= l size_x)
                     (make-array (length Y) :element-type 'cbook-symbol :initial-contents (reverse Y))
                     (if (zerop-emission-prob state (elt X l))
                         (progn (warn "dead end") (rec size_x -1 Y))
                         (let ((next_state (select-random accumA :indices-1 (list state)))
                               (Yr (select-random accumB :indices-1 (list state 0) :fixed-max +1-prob+)))
                           (if Yr
                               (rec l next_state (cons Yr Y)) ;epsilon on X
                               (let ((Yr (select-random accumB :indices-1 (list state (elt X l)))))
                                 (rec (1+ l) next_state (if (= Yr +epsilon-cbook-index+)
                                                            Y ;epsilon on Y
                                                            (cons Yr Y))))))))))
        (rec 0 init_state nil)))))
