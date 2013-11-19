;; Author: Juan Miguel Cejuela
;; Created: Wed Jul  9 18:09:49 2008 (CEST)
;; Last-Updated: 2011-08-11
;;           By: Juan Miguel Cejuela
;;     Update #: 39

(in-package :cl-hmm)

;;(declaim (optimize (speed 0) (safety 3) (compilation-speed 0) (debug 3)))
(declaim (optimize (speed 3) (safety 0)))
#+sbcl (declaim (sb-ext:muffle-conditions sb-ext:compiler-note))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun pwd (&optional (subpath ""))
    (asdf:system-relative-pathname 'cl-hmm subpath))

  (deftype cbook-state ()
    "Codebook state type. Defines in practice the maximum number of states" ;important to save memory in viterbi
    `(unsigned-byte 16))
  (deftype cbook-states ()
    "Vector of codebook states"
    `(simple-array cbook-state (*)))
  (deftype cbook-symbol ()
    "Codebook symbol type. Defines in practice the maximum number of symbols, alphabet size"
    `(unsigned-byte 8))
  (deftype cbook-alphabet ()
    "Vector of codebook symbols, alphabet codebook type"
    `(simple-array cbook-symbol (*))) ;TODO study to make it a simple-vector

  (defparameter *prob-float*
    'double-float
    "Current float precision to work with")
  (deftype prob-float ()
    "precision-probability-float, precision for probabilities"
    *prob-float*)
  (deftype prob-array (dimensions)
    `(simple-array prob-float ,dimensions))
  (defmacro prob (n)
    "Transform the number to the current float precision, *prob-float*"
    `(coerce ,n 'prob-float))

  (defconstant +0-prob+ (prob 0))
  (defconstant +1-prob+ (prob 1))
  (defconstant +most-negative-prob-float+ (symbol-value (intern (format nil "MOST-NEGATIVE-~a" *prob-float*))))
  (defconstant +very-negative-prob-float+ (/ +most-negative-prob-float+ (expt 10 8))) ;to avoid underflows
  (defconstant +most-positive-prob-float+ (symbol-value (intern (format nil "MOST-POSITIVE-~a" *prob-float*))))
  (defconstant +least-negative-prob-float+ (symbol-value (intern (format nil "LEAST-NEGATIVE-~a" *prob-float*))))
  (defconstant +least-positive-prob-float+ (symbol-value (intern (format nil "LEAST-POSITIVE-~a" *prob-float*))))

  ;;; Empty emission epsilon symbol, ε
  (defconstant +epsilon+ "")
  (defconstant +epsilon-cbook-index+ 0)

  ;;other
  (defconstant +buffer-stream-size+ 1024) ;buffer size for output random generated seqs. see hmm-run
  (defconstant +stream-correctness-size+ 4096) ;;see hmm-correctp

  ;;list of specialized hmm's classes.
  (defparameter *hmm-types* nil))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; metaclass HMMs to cover all hmm methods
(defclass hmm () nil)

;;; Define all sub HMMs through this macro: defines its multiple slot accessor, and add the name type to the global list
;;; of HMMs, *hmm-types*.
;;;
;;; multiple-accessor explanation: see with-typed-slot-values (jmc.cl.utils). Compulsory for efficiency reasons
;;;
(defmacro def-hmm-type (name direct-superclasses slot-types-name multi-accessor-name direct-slots &rest options)
  "Define a HMM class type:
  name: name of the HMM class
  direct-superclasses: superclasses, same as in DEFCLASS
  slot-types-name: name for the proper list that will store the types for the HMM class, see with-typed-slot-values (jmc.cl.utils)
  multi-accessor-name: name for the multi accessor macro, see with-typed-slot-values (jmc.cl.utils)
  direct-slots: slot definitions, same as in DEFCLASS
  options: options, same as in DEFCLASS"
  `(eval-when (:compile-toplevel :load-toplevel :execute)
     (progn
       (defclass ,name ,direct-superclasses
         ,direct-slots ,@options)
       ,@(when slot-types-name
               `((defparameter ,slot-types-name nil)
                 (funcall #'(lambda ()
                              (let ((type))
                                (dolist (slot ',direct-slots)
                                  (setf type (member :type slot))
                                  (when type (push (cons (car slot) (second type)) ,slot-types-name)))
                                (setf ,slot-types-name (nreverse ,slot-types-name)))))
                 (defmacro ,multi-accessor-name (slots instance &body body)
                   `(with-typed-slot-values ,slots ,instance ,,slot-types-name ,@body))))
       (push ',name *hmm-types*))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;
;;;; Rules to define a new hmm type
;;;;   Each new hmm or group of hmms must be declared in a new file called after it
;;;;   All hmm must be defined through the macro def-hmm-type
;;;;   Order of specification: definition, initialization methods, common methods, specific methods
;;;;


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Initialization
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; initializate-instance specification
;; (defun make-hmm(-*) ()
;; (defun make-random-hmm(-*) ()
;; (defun make-uniform-hmm(-*) ()


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Common methods
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric hmm-copy (hmm)
  (:documentation "Copy the hmm
  hmm: hmm model
  value: hmm model"))

(defgeneric hmm-correctp (hmm)
  (:documentation "T if the hmm is a correct hmm for its type, NIL otherwise
  hmm: hmm model
  value1: boolean
  value2: string / performed tests results"))

(define-condition hmm-incorrect (error)
  ((text :initarg :text :reader text)))

(defgeneric hmm-incorrect-signal (hmm))

(defgeneric hmm-compatiblep (hmm1 hmm2)
  (:documentation "T if the hmms could transform into the other because share the same characteristics. If NIL prints the reasons
  hmm: hmm model
  value: boolean"))

(defgeneric hmm-run (hmm &optional max-length)
  (:documentation "Run the hmm, ie, generate a random sequence
  hmm: hmm model
  max-length: integer / max length for the generated sequence
  verbose: boolean / if T outputs the state-pathway"))

(defgeneric hmm-no-transitions (hmm)
  (:documentation "Number of transitions
  hmm: hmm model
  value: integer"))

(defgeneric hmm-complexity (hmm)
  (:documentation "Model complexity relative to the hmm model type
  hmm: hmm model
  value: integer"))

(defgeneric cbook-encode (hmm observation)
  (:documentation "cbook-encode observation sequence of symbols to sequence of cbook-indexes"))

(defgeneric cbook-decode (hmm observation)
  (:documentation "cbook-decode observation sequence of cbook-indexes to sequence of symbols"))

(defgeneric !hmm-noisify (hmm noise)
  (:documentation "(DESTRUCTIVE) Randomly add noise to the model parameters
  @param noise: float in [0, 1]; 0 to not change the model, 1 for completely random model"))

(defgeneric hmm-save (hmm filename &optional model-spec)
  (:documentation
   "Save hmm to file

    @param hmm: hmm to save
    @param filename: filename to save the hmm to
    @param model-spec (optional): [complete|relevant] specification-way to save the model

    @return nil"))

;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Specific methods
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;...
