;; Author: Juan Miguel Cejuela
;; Created: Wed Jul  9 19:05:15 2008 (CEST)
;; Version:
;; Last-Updated: 2011-08-11
;;           By: Juan Miguel Cejuela
;;     Update #: 85
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;; Description:
;;
;; This file defines the basic HMM type, the hmm-simple class. This is divided
;; into the infinite and finite classes. They are almost the same and only
;; differ in the last steps of viterbi/forward/backward. The model
;; characteristics are well docummented in the main class definition.
;;
;; Every state keeps a bit-properties vector, namely:
;; * begin-state (true if its init prob is non zero)
;; * end
;; * silent
;; * invalid (no transitions)
;; * tied emissions
;; * fixed (the probabilities don't change, though not implemented by now)
;;
;; Every state can be identified by either an atom or a list containing an atom
;; and an index. In the latter, the atom represents the name of the group that
;; the state belongs to. Probability emissions are tied automatically by their
;; name or group-name. At the same time a state can has a label identifier that
;; represents its meaning of the context. Here is the structure of a state:
;;
;; (bit-vector (or nil atom (list atom integer)) (or nil or character))
;; (properties   state-name (group-name&index)           label
;;
;; The probabilities are saved in arrays. In order to save calculations when the
;; model is linear, two lists expressing the transition connections between the
;; states are kept. These are, the connections from state x, and the connections
;; to state x (used this in viterbi and forward).
;;

(in-package :net.ashrentum.cl-hmm)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(eval-when (:compile-toplevel :load-toplevel :execute)
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; states
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; mapping: begin | end | silent | invalid  ||  tied | fixed | x | x
  (defconstant +begin-state-flag+ 0)
  (defconstant +end-state-flag+ 1)
  (defconstant +silent-state-flag+ 2)
  (defconstant +invalid-state-flag+ 3)
  (defconstant +tied-state-flag+ 4)
  (defconstant +fixed-state-flag+ 5)
  (deftype state-properties ()
    "Bit vector for the state properties"
    `(bit-vector 8))
  (deftype state-group ()
    "Group name for a state. Group of tied emissions"
    `atom)
  (defun state-name-p (a)
    (or (typep a 'state-group) (and (listp a) (typep a 'state-group) (typep (second a) 'cbook-state))))
  (deftype state-name ()
    "Name type a state. Group name + numerical identifier"
    `(or null (satisfies state-name-p)))
  (deftype state-label ()
    "Label name for a state. Meaning/representation of the state (eg: non-coding state: N)"
    `character)

  (defconstant +label-wildcard+ #\*) ;can match any state
  (defconstant +label-null+ #\?) ;not specified

  (defun state-p (a)
    (and (listp a)
         (typep (first a) 'state-properties)
         (typep (second a) 'state-name)
         (typep (third a) 'state-label)))

  (deftype state () ;;****
    "State constituent"
    `(satisfies state-p))

  (deftype vector-states ()
    "Vector of states"
    `(simple-array state (*)))

  (defun state-properties (state)
    "Properties of the state"
    (declare (inline state-properties))
    (first state))

  (defun state-name (state)
    "Name of the state"
    (declare (inline state-name))
    (second state))

  (defun print-pretty-state-name (state-name)
    "Given the name of a state (the group and the index) concatenate that information"
    (cond
      ((null state-name) nil)
      ((listp state-name)
       (if (symbolp (car state-name))
           (concatenate 'string (symbol-name (car state-name)) (write-to-string (second state-name)))
           (concatenate 'string (write-to-string (car state-name)) (write-to-string (second state-name)))))
      (t (if (symbolp state-name)
             (symbol-name state-name)
             (write-to-string state-name)))))

  (defun state-group (state)
    "Give state's group name, if any"
    (declare (inline state-group))
    (let ((sn (state-name state)))
      (cond
        ((null sn) nil)
        ((atom sn) sn)
        (t (first sn)))))

  (deftype vector-state-groups ()
    "Vector of state-groups"
    `(simple-array state-group (*)))

  (defun state-group-index (state)
    "Group index of the state"
    (declare (inline state-group-index))
    (second (state-name state)))

  (defun state-label (state)
    "Label of the state"
    (declare (inline state-label) (optimize (speed 3) (safety 0)))
    (third state))

  (defun state-label= (label1 label2)
    "Compare two labels"
    (declare (inline state-label=) (state-label label1 label2) (optimize (speed 3) (safety 0)))
    (char= label1 label2))

  (defmacro compatible-label (state-label symbol-label)
    "Resolve label compatibility with the one given in symbol-label"
    `(or (state-label= (the state-label +label-wildcard+) (the state-label ,symbol-label))
         (state-label= (the state-label ,state-label) (the state-label ,symbol-label))))

  ;;values no-groups vector-of-groups-names vector-of-state-groups (index coded)
  (defun define-groups (states)
    (declare (vector-states states))
    (let* ((N (length states))
           (group)
           (groups (make-array N :adjustable t :fill-pointer 0))
           (grouped (make-array N :element-type 'cbook-state :initial-element 0))
           (no-groups 0)
           (pos 0)
           (nil-group))
      (declare (cbook-state N no-groups))
      (dotimes (i N (values no-groups (coerce groups '(simple-array state-group (*))) grouped))
        (setf group (state-group (aref states i)))
        (unless (and group (setq pos (position group groups)))
          (incf no-groups)
          (if group
              (vector-push group groups)
              (unless nil-group
                (setq nil-group pos)))
          (setq pos (1- no-groups)))

        (setq pos (if (and group nil-group (>= nil-group pos))
                      (1+ pos)
                      pos))
        (setq pos (coerce pos 'cbook-state))
        (setf (aref grouped i) pos))))

;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Alphabet
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  (deftype emission-symbol ()
    "Emission symbol type. Accepts any type."
    t)

  (deftype alphabet ()
    `(simple-array));; emission-symbol)) let it to be specialized

;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Model properties
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  (deftype PE-vec ()
    `(prob-array (*)))
  (deftype A-array ()
    `(prob-array (* *)))
  (deftype B-array ()
    `(prob-array (* *)))

;;; info transition lists. Vector of lists which store the index of the connected states
  (defun cbook-state-list-p (list)
    (and (listp list)
         (or (null list)
             (every #'(lambda (x) (typep x 'cbook-state)) list))))
  (deftype itrans-list ()
    "List of connected states"
    `(satisfies cbook-state-list-p))
  (deftype itrans ()
    "Vector of list of connected states"
    `(simple-array itrans-list (*)))

  (defmacro dolist-itrans ((var list &optional result) &body body)
    "Transverse a itrans-list"
    (let ((l (gensym)))
      `(do ((,l ,list (cdr ,l))
            (,var 0))
           ((null ,l) ,result)
         (declare (cbook-state ,var) (itrans-list ,l))
         (setq ,var (car ,l))
         ,@body)))

  (defun trans-array->itrans (A &optional (transpose nil))
    ;;notice that when A is transposed, the result is the info of the states that connect to
    (let* ((A. (if transpose (transpose A) A))
           (N (array-dimension A. 0))
           (itrans (make-typed-array N 'itrans-list nil))
           (l nil))
      (declare (itrans itrans))
      (dotimes (i N (coerce itrans 'itrans))
        (setq l nil)
        (dotimes (j (array-dimension A. 1))
          (unless (zerop (aref A. i j))
            (push j l)))
        (setf (aref itrans i) (nreverse l)))))

;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Parsing the model
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  (defconstant +connect-to-me+ :@) ;shortcut to specify a loop
  (defconstant +connect-to-next+ :>) ;shortcut to specify a conexion to the next state

);;;;end eval-when

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Parsing the model
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmacro sis->contribute-init (sis state-index PE) ;sis, after state init specification
  `(setf (aref ,PE ,state-index) (+ (aref ,PE ,state-index) (prob ,sis))))

;;(0 0.7 :biased 0.3) -> (0 0.7 1 0.3) and check of construction errors
(defun sas-check&fix (sas state-index S-hash) ;sas, after state transition (A) specification
  (let ((N (hash-table-count S-hash))
        (sas.))
    (if (numberp sas) ;next state very shortcut
        (list (1+ state-index) sas)
        (do ((s (first sas) (first l))
             (prob (second sas) (second l))
             (l (cddr sas) (cddr l)))
            ((null s) (nreverse sas.))
          (if (numberp s)
              (unless (and (>= s 0) (< s N)) (error "Invalid index for a state: ~a (must be 0 < x < ~a)" s N))
              (cond
                ((equal s +connect-to-me+) (setf s state-index)) ;circular shortcut
                ((equal s +connect-to-next+) (setf s (1+ state-index))) ;next state shortcut
                (t (let ((v (gethash s S-hash)))
                     (unless v (error "Not such a name for any state: ~a ~a" s v))
                     (setf s v)))))
          (push s sas.)
          (check-type prob number)
          (push prob sas.)))))

(defmacro sas->contribute-trans (sas state-index A)
  (let ((s (gensym)) (p (gensym)) (l (gensym)))
  `(do ((,s (first ,sas) (first ,l))
        (,p (second ,sas) (second ,l))
        (,l (cddr ,sas) (cddr ,l)))
       ((null ,p) ,A)
     (setf (aref ,A ,state-index ,s) (+ (aref ,A ,state-index ,s) (prob ,p))))))

(defun ses-check&fix (ses groups-emis-hash)
  (let ((emis))
    (cond
      ((null ses) nil)
      ((or (atom ses) (and (leng=2 ses) (setq emis (gethash ses groups-emis-hash))))
         (unless emis
           (setq emis (gethash ses groups-emis-hash)))
       emis)
      (t ses))))

(defmacro ses->contribute-emis (ses state-index B) ;ses, after state emission specification
  (let ((s (gensym)) (p (gensym)) (l (gensym)))
  `(do ((,s 0 (1+ ,s))
        (,p (first ,ses) (first ,l))
        (,l (cdr ,ses) (cdr ,l)))
       ((null ,p) ,B)
     (setf (aref ,B ,state-index ,s) (+ (aref ,B ,state-index ,s) (prob ,p))))))

(defmacro sms-contribute-model (state-index PE A B sis sas ses S-hash groups-emis-hash) ;sms, after state model spec
  `(progn
     (when ,sis (sis->contribute-init ,sis ,state-index ,PE))
     (sas->contribute-trans (sas-check&fix ,sas ,state-index ,S-hash) ,state-index ,A)
     (ses->contribute-emis (ses-check&fix ,ses ,groups-emis-hash) ,state-index ,B)))


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Pretty printing model properties
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun init-stats (PE S)
  (declare (PE-vec PE) (vector-states S))
  (let ((N (array-dimension  PE 0)) (init) (counter 0))
    (dotimes (i N (values counter (nreverse init)))
      (unless (zerop (aref PE i))
        (incf counter)
        (push (print-pretty-state-name (state-name (aref S i))) init)
        (push (aref PE i) init)))))

(defun trans-stats (A S &optional (transpose nil) (real-name nil))
  (declare (A-array A) (vector-states S))
  (let* ((N (array-dimension A 0)) (out) (indv-out) (A. (if transpose (transpose A) A)) (counter 0) (name))
    (declare (fixnum N) (A-array A.) (fixnum counter))
    (dotimes (i N (values (reduce #'+ out :key #'second) (nreverse out)))
      (setf indv-out nil
            counter 0)
      (dotimes (j N)
        (unless (zerop (aref A. i j))
          (incf counter)
          (setf name (if real-name (state-name (aref S j)) (print-pretty-state-name (state-name (aref S j)))))
          (unless name (setf name j))
          (push name indv-out)
          (push (aref A. i j) indv-out)))
      (setf name (if real-name (state-name (aref S i)) (print-pretty-state-name (state-name (aref S i)))))
      (unless name (setf name i))
      (push (cons name (cons counter (nreverse indv-out))) out))))

(defun emis-stats (B alphabet S &optional (real-name nil))
  (declare (A-array B) (vector-states S))
  (let ((N (array-dimension B 0))
        (M (array-dimension B 1))
        (out nil) (indv-out nil)
        (counter 0)
        (symbol)
        (name)
        (group)
        (emissions)
        (groups-emissions (make-hash-table :test 'equalp)))
    (declare (fixnum N M) (fixnum counter))
    (dotimes (i N (nreverse out))
      (setq group (state-group (aref S i)))
      (if (and (tied-state-3p i S) (setq emissions (gethash group groups-emissions)))
          (push (list (if real-name
                          (state-name (aref S i))
                          (print-pretty-state-name (state-name (aref S i))))
                      (second emissions) (first emissions)) out)
          (progn
            (setq indv-out nil
                  counter 0)
            (dotimes (j M)
              (unless (zerop (aref B i j))
                (incf counter)
                (push (if (or (characterp (setq symbol (aref alphabet j))) (numberp symbol)) symbol j) indv-out)
                (push (aref B i j) indv-out)))
            (setf name (if real-name
                           (state-name (aref S i))
                           (print-pretty-state-name (state-name (aref S i)))))
            (unless name (setf name i))
            (setq indv-out (cons name (cons counter (nreverse indv-out))))
            (push indv-out out)
            (when (tied-state-3p i S) (setf (gethash group groups-emissions) indv-out)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Definition
;;
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def-hmm-type hmm-simple (hmm) *hmm-simple-types* hmm-simple-slots
  ((name
    :initarg :name :initform nil :accessor hmm-name)
   (N ;number of states
    :initarg :N :type cbook-state :accessor hmm-no-states)
   (M ;discrete alphabet size
    :initarg :M :type cbook-symbol :accessor hmm-no-emissions)
   (V ;individual observation symbols, alphabet
    :initarg :V :type alphabet :initform (error "Must set the alphabet") :accessor hmm-alphabet)
   (V-hash ;observation symbols -> index
    :initarg :V-hash :type hash-table :accessor hmm-alphabet-hash)
   (S ;individual states
    :initarg :S :type vector-states :initform nil :accessor hmm-states)
   (S-hash ;state symbols -> index
    :initarg :S-hash :type hash-table :accessor hmm-states-hash)
   (no-groups ;groups are the states which emission probs are tied
    :type cbook-state :initform 0 :accessor hmm-no-groups)
   (groups ;group names
    :type vector-state-groups :accessor hmm-groups)
   (state-groups
    :type cbook-states :accessor hmm-state-groups)
   (PE ;initial state distribution, called PI in Rabiner //PE is PI in Phoenician
    :initarg :PE :type PE-vec :initform (error "Must set the init probabilities") :accessor hmm-init)
   (A ;state transition probability distribution
    :initarg :A :type A-array :initform (error "Must set the transition probabilities") :accessor hmm-trans)
   (iA-from ;transitions from states to states list
    :type itrans :accessor hmm-itrans-from)
   (iA-to ;transitions to states from states list
    :type itrans :accessor hmm-itrans-to)
   (B ;observation symbol probability distribution
    :initarg :B :type B-array :initform (error "Must set the emission probabilities") :accessor hmm-emis)))


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Finite and Infinite HMMs
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def-hmm-type hmm-infinite (hmm-simple) nil nil nil)

(def-hmm-type hmm-finite (hmm-simple) nil nil nil) ;TODO, currently only infinite is fully implemented

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Initialization
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod initialize-instance :after ((hmm hmm-simple) &key)
  (hmm-simple-slots (N M V V-hash S S-hash A) hmm
    (dotimes (i M) (setf (gethash (aref V i) V-hash) i)) ;V hash
    (dotimes (i N) (setf (gethash (state-name (aref S i)) S-hash) i)) ;states hash
    (multiple-value-bind (no-groups groups state-groups) (define-groups S)
                         (setf (hmm-no-groups hmm) no-groups)
                         (setf (hmm-groups hmm) groups)
                         (setf (hmm-state-groups hmm) state-groups))
    (setf (hmm-itrans-from hmm) (trans-array->itrans A))
    (setf (hmm-itrans-to hmm) (trans-array->itrans A t))
    (hmm-state-properties-set hmm)))

(defun make-hmm-simple (no-states no-emissions alphabet model &key name (alphabet-type T) (model-spec :relevant))
  "Make an hmm-simple. Two ways to specify the model parameters as follows:
  no-states:
  no-emissions:
  alphabet: list of emissions symbols (eg, '(A C G T))
  model: depending on model-spec this list has 2 different forms
  model1: (model-spec = :complete)
    states: state names and labels
    init: initial probs
    trans: trans probs
    emis: trans probs
    example: '(((:fair #\F) (:biased #\B)) (1 0) ((.95 .05) (.15 .85)) ((1/6 ...) (1/2 1/10 ...)))

  model2: (model-spec = :relevant)
  every list represent the information of an individual state, format: state-name state-label init trans emis
    example: ((:fair 0) #\F 0.95 (:fair .95 :biased .05) (1/6 1/6 1/6 1/6 1/6 1/6))"
  (let* ((N no-states)
         (M no-emissions)
         (states)
         (V (make-array M :element-type alphabet-type))
         (V-hash)
         (S (make-array N :element-type 'state :initial-element (list (make-typed-array 8 'bit 0))))
         (S-hash)
         (PE (make-typed-array N 'prob-float +0-prob+))
         (A (make-typed-array (list N N) 'prob-float +0-prob+))
         (B (make-typed-array (list N M) 'prob-float +0-prob+))
         (groups-emis-hash (make-hash-table :test 'equal))
         (type))
;;;what is different according to the way of giving the model specification
    (ecase model-spec
      (:complete
       (setf states (first model))
       (setf S-hash (make-hash-table-with-list states :test 'equal)) ;S-hash
       (dolist-index (a (second model) i) (setf (aref PE i) (prob a))) ;init
       (dolist-index (l (third model) i) ;trans
         (dolist-index (at l j)
           (setf (aref A i j) (prob at))))
       (dolist-index (l (fourth model) i) ;emis
         (dolist-index (a l j)
           (setf (aref B i j) (prob a)))))
      (:relevant
       (setf S-hash (make-hash-table-with-list (map 'list #'car model) :test 'equal)) ;S-hash
       (let ((name) (group) (emis))
         (dolist-index (l model i)
           (setq name (first l)
                 group (if (listp name) (car name) name))
           (setq emis (gethash group groups-emis-hash))
           (unless emis
             (setf emis (fifth l)
                   (gethash group groups-emis-hash) emis))
           (push (list name (second l)) states)
           (sms-contribute-model i PE A B (third l) (fourth l) emis S-hash groups-emis-hash)))
         (setf states (nreverse states))))
;;;what is the same
    (dolist-index (symbol alphabet i) ;alphabet
      (setf (aref V i) symbol))
    (setf V-hash (make-hash-table-with-list alphabet))
    (dolist-index (e states i) ;the states
      (setf (aref S i) (cons (make-typed-array 8 'bit 0) (if (listp e) e (list e +label-null+)))))
    (!normalize-vector PE) ;normalize probs
    (!normalize-matrix A)
    (!normalize-matrix B)

    (setf type 'hmm-infinite) ;;TODO downgrade. Currently only infinite is supported
    (make-instance type :name name :N N :M M :V V :V-hash V-hash :S S :S-hash S-hash :PE PE :A A :B B)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun make-random-hmm-simple
    (no-states no-emissions &key (eccentricity 2) states alphabet name (alphabet-type T))
  "Make a hmm-simple with no biased info
  no-states
  no-emissions
  eccentricity: real / random eccentricity. The bigger, the more dispair. If 0, all the prob are equal
  states: list / optional list of states info (eg, ((:fair #\F) (:biased #\B)))
  alphabet: list / optional alphabet (eg, '(A C G T))
  name: optional name
  alphabet-type: type of the symbols"
  (make-hmm-simple no-states no-emissions
                   (if alphabet
                       alphabet
                       (do ((a) (i 0 (1+ i)))
                           ((= i no-emissions) (nreverse a))
                         (push i a)))
                   (list
                    (if states
                        states
                        (do ((a) (i 0 (1+ i)))
                            ((= i no-states) (nreverse a)) (push i a)))
                    (make-list-meval no-states (expt (random 1.0) eccentricity))
                    (make-list-meval no-states (make-list-meval no-states (expt (random 1.0) eccentricity)))
                    (make-list-meval no-states (make-list-meval no-emissions (expt (random 1.0) eccentricity))))
                   :name name :alphabet-type alphabet-type :model-spec :complete))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod print-object ((object hmm-simple) stream)
  (declare (stream stream))
  (labels ((print-A (A S)
             (let ((out ""))
               (dolist-index (trans (multiple-value-bind (a b) (trans-stats A S)
                                      (setf out (format nil " (~a): ~&" a)) b) i out)
                 (setf out
                       (concatenate 'string out (format nil "~5T~a (~a):  ~{~a:~3$~^  ~}~%"
                                                        (car trans) (second trans) (cddr trans)))))))
           (print-B (B alphabet S)
             (let ((out ""))
               (dolist-index (emis (emis-stats B alphabet S) i out)
                 (setf out
                       (concatenate 'string out (format nil "~5T~a (~a):  ~{~#[~;---> ~a~:;~a:~3$~^  ~]~}~%"
                                                        (car emis) (second emis) (cddr emis))))))))
  (print-unreadable-object (object stream :type t)
    (hmm-simple-slots (name N M V S PE A B) object
      (let ((no-begins) (begins))
        (multiple-value-setq (no-begins begins) (init-stats PE S))
        (format stream "~a (~d,~d)~%  S: ~a~2%  V: ~a~2%  PE (~a):  ~{~a:~3$~^  ~}~2%  A~a~%  B:~&~a"
                name N M S V
                no-begins
                begins
                (print-A A S)
                (print-B B V S)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Specific Procedures
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;; codebook, translation between alphabet symbols and their index
(macrolet ((operation (source)
             `(dotimes (i len output)
                (setf (aref output i) ,source)))

           (with- ((&rest slots) output-type &body body)
             `(hmm-simple-slots ,slots hmm
                (let* ((len (the fixnum (length sequence)))
                       (output (make-array len :element-type ,output-type)))
                  (declare ((simple-array) output) (fixnum len))
                  ,@body))))

  (defun cbook (hmm sequence)
    "Codebook, from an observation of symbols to their indexes"
    (declare (optimize (speed 3) (safety 0)) ((vector) sequence))
    (with- (V-hash) 'cbook-symbol
           (operation (gethash (aref sequence i) V-hash))))

  (defun cbook-indexes (hmm sequence)
    "Codebook, from an observation of the indexes to their symbols"
    (declare (optimize (speed 3) (safety 0)) (cbook-alphabet sequence))
    (with- (V) (array-element-type V)
           (operation (aref V (aref sequence i))))))

(defun cbook-list (hmm sequences)
  "Translate the list of sequences in a list of sequences index-coded"
  (labels ((aux (hmm sequences cbooks)
             (cond
               ((null sequences) (nreverse cbooks))
                (t (aux hmm (cdr sequences) (cons (cbook hmm (car sequences)) cbooks))))))
    (aux hmm sequences nil)))

(defun hmm-state-labels (hmm) ;review if make it common for all
  "Create a vector with the state's labels"
  (hmm-simple-slots (N S) hmm
    (let ((out (make-typed-array N 'state-label +label-null+)))
      (dotimes (i N out)
        (setf (aref out i) (state-label (aref S i)))))))

;;; not quite nice form
(defmacro hmm-simple-alter-model (N M PE A B alpha)
  "Alter, add noise, randomly the model
  alpha: confidence in the current model (0 to 1) 1 to don't change the model"
  `(unless (= +1-prob+ ,alpha)
     (!normalize-vector
      (!combine-float-vectors
       ,PE (!normalize-vector (make-random-vector ,N +1-prob+ 'prob-float)) ,alpha 'prob-float t))
     (!normalize-matrix
      (!combine-float-matrices
       ,A (!normalize-matrix (make-random-matrix (list ,N ,N) +1-prob+ 'prob-float)) ,alpha 'prob-float t))
     (!normalize-matrix
      (!combine-float-matrices
       ,B (!normalize-matrix (make-random-matrix (list ,N ,M) +1-prob+ 'prob-float)) ,alpha 'prob-float t))))

;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; State properties
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;; REAL - First define the methods that really clarify the states's characteristics, only used in the hmm initializations
;;;; Then interface wrappers that only check the bit-vectors of the states

(defun begin-state-realp (index PE)
  (/= +0-prob+ (aref PE index))) ;not zero in the initial probs

(defun end-state-realp (index A)
  (= +1-prob+ (aref A index index)))

(defun silent-state-realp (index B)
  (let* ((len (array-dimension B 1))
         (init (* len index))
         (top (+ init len -1)))
    (declare (A-array B))
    (loop for j from init to top do ;ie, zero for all emis
         (unless (zerop (row-major-aref B j)) (return nil))
         finally (return t))))

(defun invalid-state-realp (index A)
  (let* ((len (array-dimension A 1))
         (init (* len index))
         (top (+ init len -1)))
    (declare (A-array A))
    (loop for j from init to top do ;ie, zero for all trans
         (unless (zerop (row-major-aref A j)) (return nil))
         finally (return t))))

(defun tied-state-realp (stateindex states)
  (let ((group (state-group (aref states stateindex)))
        (tied nil)
        (groupp nil)
        (no-states (array-dimension states 0)))
    (when group
      (dotimes (index no-states (nreverse tied))
        (setq groupp (state-group (aref states index)))
        (when (and (equalp group groupp) (/= index stateindex))
          (push index tied))))))

(defun state-property-set (states stateindex flag)
  (declare (inline state-property-set))
  (setf (bit (state-properties (aref states stateindex)) flag) 1))

(defun hmm-state-properties-set (hmm) ;initialize the properties
  (macrolet ((seter (flag)
               `(state-property-set S i ,(symbol-value flag))))
    (hmm-simple-slots (N S PE A B) hmm
      (dotimes (i N)
        (when (begin-state-realp i PE) (seter +begin-state-flag+))
        (when (end-state-realp i A) (seter +end-state-flag+))
        (when (silent-state-realp i B) (seter +silent-state-flag+))
        (when (invalid-state-realp i A) (seter +invalid-state-flag+))
        (when (tied-state-realp i S) (seter +tied-state-flag+))))))

(defun get-statename-index (statename hmm)
  (declare (inline get-statename-index) (state-name statename) (hmm-simple hmm))
  "Get the index for the given state and hmm"
  (gethash statename (hmm-states-hash hmm)))

(defmacro state-predicates (name flag doc)
  "Define diverse predicators for the same property depending on the avalaible information"
  (labels ((name-def (specializer)
             (intern (concatenate 'string name "-STATE-" specializer)))
           (state-property-get (states stateindex flag)
             `(when (= 1 (bit (state-properties (aref ,states ,stateindex)) ,flag)) t))) ;otherwise nil
    (let ((gdoc doc) (gflag flag))
      `(progn
         (defun ,(name-def "P") (statename hmm) ;given state name and hmm
           (declare (state-name statename))
           ,(format nil "~a. Required state name and hmm" gdoc)
           (let ((stateindex (get-statename-index statename hmm))
                 (states (hmm-states hmm)))
             (declare (vector-states states))
             (if stateindex
                 (values ,(state-property-get 'states 'stateindex gflag) stateindex)
                 (values nil nil))))

         (defun ,(name-def "2P") (stateindex hmm) ;given state index and hmm
           ,(format nil "~a. Required state index and hmm" gdoc)
           (declare (cbook-state stateindex))
           (let ((states (hmm-states hmm)))
             (if (and (< stateindex (the cbook-state (hmm-no-states hmm))) (>= stateindex 0))
                 (values ,(state-property-get 'states 'stateindex gflag) stateindex)
                 (values nil nil))))

         (defun ,(name-def "3P") (stateindex states) ;given state index and states
           (declare (cbook-state stateindex) (vector-states states))
           ,(format nil "~a. Required state index and states" gdoc)
           (declare (cbook-state stateindex))
           ,(state-property-get 'states 'stateindex gflag))))))


(state-predicates "BEGIN" +begin-state-flag+ "Return T if state is BEGIN, and NIL otherwise")
(state-predicates "END" +end-state-flag+ "Return T if state is END, and NIL otherwise")
(state-predicates "SILENT" +silent-state-flag+ "Return T if state is SILENT, and NIL otherwise")
(state-predicates "INVALID" +invalid-state-flag+ "Return T if state is invalid (ie, no trans), and NIL otherwise")
(state-predicates "TIED" +tied-state-flag+ "Return the list of tied states if state is tied, and NIL otherwise")
(state-predicates "FIXED" +fixed-state-flag+ "Return T if state is fixed (ie, can't be changed in the estimation process), and NIL otherwise")


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Common Methods
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod hmm-copy ((hmm hmm-simple))
  (let ((N (hmm-no-states hmm))
        (M (hmm-no-emissions hmm)))
    (make-instance (type-of hmm)
                   :N N
                   :M M
                   :V (copy-seq (hmm-alphabet hmm))
                   :V-hash (make-hash-table :size M :test 'equalp)
                   :S (copy-seq (hmm-states hmm))
                   :S-hash (make-hash-table :size N :test 'equalp)
                   :PE (copy-seq (hmm-init hmm))
                   :A (copy-matrix (hmm-trans hmm))
                   :B (copy-matrix (hmm-emis hmm)))))

(defmethod hmm-correctp ((hmm hmm-simple))
  "T if all accumulative probabilities are all equal to 1, and all states are reached"
  (let ((buffer (make-array +stream-correctness-size+ :element-type 'character :adjustable t :fill-pointer 0))
        (correct t))
    (labels ((valid (accum-prob)
               (cond
                 ((consp accum-prob) (setq correct (and correct t)) t)
                 ((eql accum-prob nil) (setq correct nil) nil)
                 (t (let ((c (< 0.99994d0 accum-prob 1.00005d0))) ;because of imprecision issues
                      (setq correct (and correct c))
                      c))))
             (print-verdict (stream prefix ok)
               (if ok (format stream "~a... OK~%" prefix) (format stream "~a... WRONG~%" prefix)))
             (print-list (list)
               (dolist (item list)
                 (format buffer "~a" item))))
      (macrolet ((check-array ((index limit) prefix source)
                   `(let ((wrong))
                      (dotimes (,index ,limit)
                        (unless (valid ,source)
                          (push (print-verdict nil (format nil "     ~d" ,index) nil) wrong)))
                      (print-verdict buffer ,prefix (not wrong))
                      (print-list (nreverse wrong)))))
        (let ((N (hmm-no-states hmm)) (M (hmm-no-emissions hmm))
              (aPE (accum-array (hmm-init hmm) 1 prob-float))
              (aA (accum-array (hmm-trans hmm) 2 prob-float))
              (aB (accum-array (hmm-emis hmm) 2 prob-float))
              (iA-to (hmm-itrans-to hmm)))
          (print-verdict buffer "Init" (valid (aref aPE (1- N))))
          (check-array (i N) "Trans" (aref aA i (1- N)))
          (check-array (i N) "Emis" (aref aB i (1- M)))
          (check-array (i (length iA-to)) "All Reached" (aref iA-to i)))))
    (values correct buffer)))

(defmethod hmm-incorrect-signal ((hmm hmm-simple))
  (multiple-value-bind (a b) (hmm-correctp hmm)
    (unless a (error 'hmm-incorrect :text b))))

(defmethod hmm-compatiblep ((hmm1 hmm-simple) (hmm2 hmm-simple))
  (if (equal (type-of hmm1) (type-of hmm2))
      (if (= (slot-value hmm1 'N) (slot-value hmm2 'N))
          (if (= (slot-value hmm1 'M) (slot-value hmm2 'M))
              T
              (prog1 nil (princ "Incompatible alphabet sizes")))
          (prog1 nil (princ "Incompatible nº states")))
      (prog1 nil (princ "Incompatible type of hmms"))))

(defmethod hmm-run ((hmm hmm-simple) &optional (max-length +buffer-stream-size+))
  (declare (optimize (speed 3) (safety 0)) (fixnum max-length))
  (restart-case (hmm-incorrect-signal hmm)
    (continue-anyway () nil))
  (hmm-simple-slots (N M V S) hmm
    (let ((A-accu (accum-array (slot-value hmm 'A) 2 prob-float))
          (B-accu (accum-array (slot-value hmm 'B) 2 prob-float))
          (buffer (make-array max-length :fill-pointer 0 :element-type (array-element-type V)))
          (path (make-array max-length :fill-pointer 0 :element-type 'cbook-state))
          (path-labels (make-array max-length :fill-pointer 0 :element-type 'state-label))
          (state-labels (hmm-state-labels hmm)))
      (do* ((cur-state
             (array-search (random +1-prob+) (accum-array (slot-value hmm 'PE) 1 prob-float))
             (- (the fixnum (array-search (random +1-prob+) A-accu index-A (the fixnum (+ index-A N)))) index-A))
            (index-A
             (* cur-state N) (* cur-state N))
            (index-B
             (* cur-state M) (* cur-state M))
            (l 0))
           ((or (= l max-length) (end-state-3p cur-state S)) (values buffer path-labels path))
        (declare (fixnum cur-state index-A index-B l))
        (unless (silent-state-3p cur-state S)
          (vector-push cur-state path)
          (vector-push (aref state-labels cur-state) path-labels)
          (vector-push (aref V (- (the fixnum (array-search (random +1-prob+) B-accu index-B (the fixnum (+ index-B M))))
                                  index-B)) buffer)
          (incf l))))))

(defmethod hmm-no-transitions ((hmm hmm-simple))
  (declare (inline hmm-no-transitions))
  (hmm-simple-slots (A S) hmm
    (values (trans-stats A S))))

;;;this definition can change
(defmethod hmm-complexity ((hmm hmm-simple))
  (declare (inline hmm-complexity))
  (hmm-simple-slots (M) hmm
    (* M (hmm-no-transitions hmm))))
