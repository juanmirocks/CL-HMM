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
   (logPE ;PE in log space
     :type PE-vec :accessor hmm-init-log)
   (A ;state transition probability distribution
    :initarg :A :type A-array :initform (error "Must set the transition probabilities") :accessor hmm-trans)
   (logA ;A in log space
     :type A-array :accessor hmm-trans-log)
   (iA-from ;list, transitions from states to states
    :type itrans :accessor hmm-itrans-from)
   (iA-to ;list, transitions to states from states
    :type itrans :accessor hmm-itrans-to)
   (B ;left&right pair observation probability distribution
    :initarg :B :type B-2streams-array :initform (error "Must set the pair emission probabilities") :accessor hmm-emis)
   (logB ;B in log space
     :type B-2streams-array :accessor hmm-emis-log)))

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
  (phmm-slots (S S-hash N L L-hash L-size R R-hash R-size) hmm

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

    (reset-instance hmm)
    ;;(hmm-state-properties-set hmm) TODO
    ))

(defmethod reset-instance ((hmm phmm))
  "Use when training, to reconfigure the parameters"
  (phmm-slots (PE A B) hmm
    (setf (hmm-init-log hmm) (log-array PE))
    (setf (hmm-trans-log hmm) (log-array A))
    (setf (hmm-emis-log hmm) (log-array B))
    (setf (hmm-itrans-from hmm) (trans-array->itrans A))
    (setf (hmm-itrans-to hmm) (trans-array->itrans A t))))

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

  (when (or (find 0 L-list) (find 0 R-list)) (error "The alphabets cannot contain the special 0 symbol reserved for epsilon"))
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
    (N L-size R-size &key (eccentricity 2) (states (range N)) (L-list (range (1+ L-size) :min 1)) (R-list (range (1+ R-size) :min 1)) name (L-alphabet-type T) (R-alphabet-type T))
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
    (N L-size R-size &key (states (range N)) (L-list (range (1+ L-size) :min 1)) (R-list (range (1+ R-size) :min 1)) name (L-alphabet-type T) (R-alphabet-type T))
  "Same as make-random-hmm-simple but with an uniform distribution for all probabilities, eccentricity 0"
  (make-random-phmm N L-size R-size :eccentricity 0 :states states :L-list L-list :R-list R-list :name name :L-alphabet-type L-alphabet-type :R-alphabet-type R-alphabet-type))


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Common Methods
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod hmm-correctp ((hmm phmm))
  (let ((buffer (make-array +stream-correctness-size+ :element-type 'character :adjustable t :fill-pointer 0))
        (correct t))

    (labels ((existsNegative (matrix)
               (dotimes (i (array-total-size matrix) nil)
                 (if (< (row-major-aref matrix i) 0) (return t))))
             (=p (expected actual &optional (accepted-error-range 1d-5))
               (< (- expected accepted-error-range) actual (+ expected accepted-error-range))))
      (macrolet ((assume (condition &rest errorMsgArguments)
                   `(unless ,condition (setq correct nil)
                            (format buffer ,@(cons (concatenate 'string (car errorMsgArguments) "~%") (cdr errorMsgArguments))))))

        (phmm-slots (N PE A B) hmm
          (assume (not (existsNegative PE)) "Negative value in PE")
          (assume (=p 1 (areduce #'+ PE)) "PE is not a valid probability distribution, total value: ~f" (areduce #'+ PE))
          (assume (not (existsNegative A)) "Negative value in A")
          (assume (=p N (areduce #'+ A)) "A is not a valid probability distribution, total value: ~f" (areduce #'+ A))
          (assume (not (existsNegative B)) "Negative value in B")
          (assume (=p N (areduce #'+ B)) "B is not a valid probability distribution, total value: ~f" (areduce #'+ B))
          (loop for i below N do (assume (zerop (aref B i 0 0)) "Non-zero! b_~d(epsilon, epsilon) = ~f" i (aref B i 0 0)))))

      (values correct buffer))))

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
  (declare (optimize (speed 3)) (simple-vector x))
  (let* ((size_x (the fixnum (length x)))
         (out_x (make-svector size_x)))
    (phmm-slots (L-hash) hmm
      (dotimes (i size_x out_x)
        (setf (svref out_x i) (gethash (svref x i) L-hash))))))

(defmethod cbook-encode-right ((hmm phmm) y)
  "cbook-encode right string, that is, the output of a pair observation"
  (declare (optimize (speed 3)) (simple-vector y))
  (let* ((size_y (the fixnum (length y)))
         (out_y (make-svector size_y)))
    (phmm-slots (R-hash) hmm
      (dotimes (i size_y out_y)
        (setf (svref out_y i) (gethash (aref y i) R-hash))))))

(defmethod cbook-encode ((hmm phmm) observation)
  "cbook-encode pair observation"
  (declare (optimize (speed 3)))
  (let* ((x (first observation))
         (y (second observation)))
    (list (cbook-encode-left hmm x) (cbook-encode-right hmm y))))

(defmethod cbook-decode-left ((hmm phmm) x)
  "cbook-decode left string"
  (declare (optimize (speed 3) (safety 3)) (simple-vector x))
  (let* ((size_x (length x)))
    (phmm-slots (L) hmm
      (let ((out_x (make-array size_x :element-type (array-element-type L) :fill-pointer nil)))
         (dotimes (i size_x out_x)
           (setf (aref out_x i) (aref L (1- (svref x i))))))))) ;-1 since in L&R epsilon is not accounted for

(defmethod cbook-decode-right ((hmm phmm) y)
  "cbook-decode right string"
  (declare (optimize (speed 3) (safety 3)) (simple-vector y))
  (let* ((size_y (length y)))
    (phmm-slots (R) hmm
      (let ((out_y (make-array size_y :element-type (array-element-type R) :fill-pointer nil)))
         (dotimes (i size_y out_y)
           (setf (aref out_y i) (aref R (1- (svref y i))))))))) ;-1 since in L&R epsilon is not accounted for

(defmethod cbook-decode ((hmm phmm) observation)
  "cbook-decode pair observation"
  (declare (optimize (speed 3) (safety 0)))
  (let* ((x (first observation))
         (y (second observation)))
    (list (cbook-decode-left hmm x) (cbook-decode-right hmm y))))

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

(defmacro semiring (space &body body)
  `(let ((ZERO    (if (eq ,space :log) +LOGZERO+      +0-prob+))
         ;;TODO define incf-like setter with SUM (it's both cleaner and faster for array access)
         (SUM     (if (eq ,space :log) 'log+          '+))
         (MINUS   (if (eq ,space :log) 'log-          '-))
         (SUMF    (if (eq ,space :log) 'logincf       'incf))
         ;;TODO log* could be replaced perhaps for simple + and just check underflow-sensitive places
         (MUL     (if (eq ,space :log) 'log*          '*))
         (ONE     (if (eq ,space :log) +0-prob+       +1-prob+))
         ;;TODO this is more than a semiring, change the name
         (DIV     (if (eq ,space :log) '-             '/)))
     `(when ,ONE (,DIV ,ONE ,ONE)) ;ignore unused
     ,@body))

(defmacro alpha[] (space dim1 dim2 dim3)
  "Alpha matrix accessor. The matrix is assumed to be called exactly 'alpha'"
  `(if (or (< ,dim2 0) (< ,dim3 0))
       ,(if (eq space :log) +LOGZERO+ +0-prob+)
       (aref alpha ,dim1 ,dim2 ,dim3)))

(defmacro cbref1 (seq i)
  "1-indexed cbook-encoded input sequence. If 0, return epsilon's index"
  (declare (optimize (speed 3) (safety 0)) (simple-vector seq) (fixnum i))
  `(if (zerop ,i)
       +epsilon-cbook-index+
       (svref ,seq (1- ,i))))

(defmacro define-forward (space)
  (semiring space
    `(defmethod ,(if (eq space :log) 'forward-log 'forward) ((hmm phmm) obs-c)
       ,(format nil "
Forward algorithm in ~a space.

@param hmm: pair hidden markov model
@param obs-c: cbook-encoded pair observation, list of 2 elements

@return (1) ~aprobability of observation pair
@return (2) alpha 3d matrix
" space (if (eq space :log) "log-" ""))

       (declare (optimize (speed 3) (safety 0)))
       (phmm-slots (N iA-to) hmm
         (let* ((PE (the PE-vec           (slot-value hmm ,(if (eq space :log) ''logPE ''PE))))
                (A  (the A-array          (slot-value hmm ,(if (eq space :log) ''logA  ''A))))
                (B  (the B-2streams-array (slot-value hmm ,(if (eq space :log) ''logB  ''B))))
                (x (first obs-c))
                (y (second obs-c))
                (size_x (length x))
                (size_y (length y))
                (alpha (the (prob-array (* * *)) (make-typed-array `(,N ,(1+ size_x) ,(1+ size_y)) 'prob-float ,ZERO))))
           (declare (simple-vector x y) (fixnum size_x size_y))

           ;; Initialization
           ;; -------------------------------------------------------------------------
           (loop for j below N do
                (when (> size_x 0) (setf (aref alpha j 1 0) (,MUL (aref PE j) (aref B j (svref x 0) +epsilon-cbook-index+))))
                (when (> size_y 0) (setf (aref alpha j 0 1) (,MUL (aref PE j) (aref B j +epsilon-cbook-index+ (svref y 0)))))
                (when (and (> size_x 0) (> size_y 0))
                  (loop for i in (aref iA-to j)
                     with diag :of-type prob-float = (aref PE j)
                     with l-1  :of-type prob-float = ,ZERO
                     with r-1  :of-type prob-float = ,ZERO
                     do
                       (setf l-1 (,SUM l-1 (,MUL (aref A i j) (aref alpha i 0 1))))
                       (setf r-1 (,SUM r-1 (,MUL (aref A i j) (aref alpha i 1 0))))
                     finally
                       (setf (aref alpha j 1 1)
                             (,SUM (,SUM
                                    (,MUL diag (aref B j (svref x 0) (svref y 0)))
                                    (,MUL l-1  (aref B j (svref x 0) +epsilon-cbook-index+)))
                                   (, MUL r-1  (aref B j +epsilon-cbook-index+ (svref y 0))))))))

           ;; Induction
           ;; -------------------------------------------------------------------------
           (loop for l from 0 to size_x do
                (loop for r from 0 to size_y do
                     (when (<= 2 (max l r))
                       (loop for j below N do
                          ;;in contrast with the semicode, traverse iA-to at the upper level and only once since this operation is costly
                            (loop for i in (aref iA-to j)
                               with diag :of-type prob-float = ,ZERO
                               with l-1  :of-type prob-float = ,ZERO
                               with r-1  :of-type prob-float = ,ZERO
                               do
                                 (setf diag (,SUM diag (,MUL (aref A i j) (alpha[] ,space i (1- l) (1- r)))))
                                 (setf l-1  (,SUM l-1  (,MUL (aref A i j) (alpha[] ,space i (1- l) r    ))))
                                 (setf r-1  (,SUM r-1  (,MUL (aref A i j) (alpha[] ,space i l      (1- r)))))

                               finally
                                 (setf (aref alpha j l r)
                                       (,SUM (,SUM
                                              (,MUL diag (aref B j (cbref1 x l)          (cbref1 y r)))
                                              (,MUL l-1  (aref B j (cbref1 x l)          +epsilon-cbook-index+)))
                                             (, MUL  r-1  (aref B j +epsilon-cbook-index+ (cbref1 y r))))))))))

           ;; Termination
           ;; -------------------------------------------------------------------------
           (values
            (loop for j below N with ret :of-type prob-float = ,ZERO do (setf ret (,SUM ret (aref alpha j size_x size_y))) finally (return ret))
            (the (prob-array (* * *)) alpha)))))))

(define-forward :log)
(define-forward :probability)

(defmacro beta[] (space dim1 dim2 dim3)
  "Beta matrix accessor. The matrix is assumed to be called exactly 'beta'"
  `(if (or (> ,dim2 size_x) (> ,dim3 size_y))
       ,(if (eq space :log) +LOGZERO+ +0-prob+)
       (aref beta ,dim1 ,dim2 ,dim3)))

(defmacro cbref1-beta (seq i)
  "1-indexed cbook-encoded input sequence. if i >= length(seq), return epsilon's index.
   Note: we don't index by (1-) since the function is already called (in backward) with the index - 1 (for efficiency)"
  `(if (= ,i (length ,seq))
       +epsilon-cbook-index+
       (svref ,seq ,i)))

(defmacro define-backward (space)
  (semiring space
    `(defmethod ,(if (eq space :log) 'backward-log 'backward) ((hmm phmm) obs-c)
       ,(format nil "
Backward algorithm in ~a space.

@param hmm: pair hidden markov model
@param obs-c: cbook-encoded pair observation, list of 2 elements

@return generated beta 3d matrix
" space)

       (declare (optimize (speed 3) (safety 0)))
       (phmm-slots (N iA-from) hmm
         (let* ((A  (the A-array          (slot-value hmm ,(if (eq space :log) ''logA  ''A))))
                (B  (the B-2streams-array (slot-value hmm ,(if (eq space :log) ''logB  ''B))))
                (x (first obs-c))
                (y (second obs-c))
                (size_x (length x))
                (size_y (length y))
                (beta (make-typed-array `(,N ,(1+ size_x) ,(1+ size_y)) 'prob-float ,ZERO)))
           (declare (simple-vector x y) (fixnum size_x size_y))

           ;;Initialization
           ;; -------------------------------------------------------------------------
           (loop for i below N do
                (setf (aref beta i size_x size_y) ,ONE))

           ;;Induction
           ;; -------------------------------------------------------------------------
           (loop for l from size_x downto 0 do
                (loop for r from size_y downto 0 do
                     (when (not (and (= l size_x) (= r size_y)))
                       (loop for i below N do
                            (setf (aref beta i l r)
                                  (loop for j in (aref iA-from i)
                                     with accum :of-type prob-float = ,ZERO
                                     do
                                       (setf accum
                                             (,SUM accum
                                                   (,MUL (aref A i j)
                                                         (,SUM (,SUM
                                                                (,MUL (beta[] ,space j (1+ l) (1+ r)) (aref B j (cbref1-beta x l)     (cbref1-beta y r)))
                                                                (,MUL (beta[] ,space j (1+ l) r     ) (aref B j (cbref1-beta x l)     +epsilon-cbook-index+)))
                                                               ( ,MUL (beta[] ,space j l      (1+ r)) (aref B j +epsilon-cbook-index+ (cbref1-beta y r)))))))
                                     finally (return accum)))))))

           (the (prob-array (* * *)) beta))))))

(define-backward :log)
(define-backward :probability)

;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Run & Translations
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod hmm-run ((phmm phmm) &optional max-length)
  (if (not max-length) (error "max-length must be given"))

  (phmm-slots (N L-size R-size PE A B) phmm
    (let* ((init_state (select-random (accum-array PE 1 prob-float)))
           (accumA (accum-array A 2 prob-float))
           (accumB (accum-array B 3 prob-float))
           (accumBL (let ((accumBLtemp (make-typed-array (list N (1+ L-size)) 'prob-float +0-prob+)))
                      (dotimes (in N)
                        (dotimes (il (1+ L-size))
                          (setf (aref accumBLtemp in il) (aref accumB in il R-size))))
                      (accum-array accumBLtemp 2 prob-float)))
           (length (1+ (random max-length))))
      (labels ((end-state-p (state) (or (= 1 (aref A state state)) (not (select-random accumA :indices-1 (list state)))))
               (rec (state size_x x size_y y)
                 (declare (fixnum size_x size_y) (list x y))
                 (if (or (end-state-p state) (and (or (= length size_x) (= length size_y)) (not (zerop size_x)) (not (zerop size_y))))
                     (list
                      (make-array size_x :fill-pointer nil :initial-contents (reverse x))
                      (make-array size_y :fill-pointer nil :initial-contents (reverse y)))
                     (let* ((next_state (select-random accumA :indices-1 (list state)))
                            (xl (select-random accumBL :indices-1 (list state)))
                            (xl_not_epsilon (/= xl +epsilon-cbook-index+))
                            (yr (the fixnum (select-random accumB :indices-1 (list state xl))))
                            (yr_not_epsilon (/= yr +epsilon-cbook-index+)))
                       (cond
                         ((and xl_not_epsilon yr_not_epsilon) (rec next_state (1+ size_x) (cons xl x) (1+ size_y) (cons yr y)))
                         (xl_not_epsilon (rec next_state (1+ size_x) (cons xl x) size_y y))
                         (yr_not_epsilon (rec next_state size_x x (1+ size_y) (cons yr y)))
                         (t (error "Impossible case, epsilon|epsilon")))))))

        (rec init_state 0 nil 0 nil)))))

(defmethod hmm-translate ((phmm phmm) X)
  "Translate X input sequence (left stream) -> Y output sequence (right stream)

  The Y translation is non-deterministic and is constructed following a random path
  determined by the transition and emission distributions of the model. Some randomly-generated
  paths may lead to dead ends, that is, where the X is not completely generated.

  @param phmm: model with probability distributions to use for translation
  @param X: cbook-encoded left-stream input sequence

@return(0) Y: a translation of X, cbook-encoded right-stream output sequence"

  (declare (optimize (speed 3) (safety 0)) (simple-vector X))

  (phmm-slots (PE A B) phmm
    (let ((size_x (length X))
          (init_state (select-random (accum-array PE 1 prob-float)))
          (accumA (accum-array A 2 prob-float))
          (accumB (accum-array B 3 prob-float)))
      (labels ((zerop-emission-prob (state l) (not (select-random accumB :indices-1 (list state l))))
               (rec (l state Y)
                 (declare (list Y))
                 (if (= l size_x)
                     (make-array (length Y) :fill-pointer nil :initial-contents (reverse Y))
                     (if (zerop-emission-prob state (svref X l))
                         (progn (warn "dead end") (rec size_x -1 Y))
                         (let ((next_state (select-random accumA :indices-1 (list state)))
                               (Yr (the fixnum (select-random accumB :indices-1 (list state +epsilon-cbook-index+) :fixed-max +1-prob+))))
                           (if (>= Yr 0)
                               (rec l next_state (cons Yr Y)) ;epsilon on X
                               (let ((Yr (the fixnum (select-random accumB :indices-1 (list state (svref X l))))))
                                 (rec (1+ l) next_state (if (= Yr +epsilon-cbook-index+)
                                                            Y ;epsilon on Y
                                                            (cons Yr Y))))))))))
        (rec 0 init_state nil)))))

(defmethod hmm-translate-viterbi ((phmm phmm) X &optional (hmm-left (hmm-left phmm)))
  (let ((viterbi-path (viterbi-log hmm-left X)))
    (phmm-slots (B) phmm
      (let ((size_x (length X))
            (accumB (accum-array B 3 prob-float)))
        (labels ((zerop-emission-prob (state l) (not (select-random accumB :indices-1 (list state l))))
                 (rec (l states-path Y)
                   (declare (list Y))
                   (if (= l size_x)
                       (make-array (length Y) :fill-pointer nil :initial-contents (reverse Y))
                       (let ((state (aref states-path l)))
                         (if (zerop-emission-prob state (svref X l))
                             (progn (warn "dead end") (rec size_x -1 Y))
                             (let ((Yr (select-random accumB :indices-1 (list state +epsilon-cbook-index+) :fixed-max +1-prob+)))
                               (if (>= Yr 0)
                                   ;;As it is now, with epsilon on X, the transition is forced to
                                   ;;state in the same state which may be illegal if P(aii) == 0
                                   (rec l states-path (cons Yr Y)) ;epsilon on X
                                   (let ((Yr (select-random accumB :indices-1 (list state (svref X l)))))
                                     (rec (1+ l) states-path (if (= Yr +epsilon-cbook-index+)
                                                                 Y ;epsilon on Y
                                                                 (cons Yr Y)))))))))))
          (rec 0 viterbi-path nil))))))

(defmethod hmm-left ((phmm phmm))
  (phmm-slots (N S L L-size R-size name PE A B) phmm
    (let ((states (loop for i below N collect (list (state-name (aref S i)) (state-label (aref S i)))))
          (newB (make-typed-array (list N L-size) 'prob-float +0-prob+)))
      (loop for i below N do
           (loop for l from 1 to L-size do
                (setf (aref newB i (1- l)) (loop for r from 0 to R-size sum (aref B i l r)))))
      (make-hmm-simple N L-size (array->list L)
                       `(,states
                         ,(array->list PE)
                         ,(array->list A)
                         ,(array->list newB))
                       :name name :model-spec :complete))))
