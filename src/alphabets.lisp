;; Author: Juan Miguel Cejuela
;; Created: Wed Jul  9 19:10:30 2008 (CEST)
;; Last-Updated: 2011-08-11
;;           By: Juan Miguel Cejuela
;;     Update #: 9
;; Keywords:

(in-package :cl-hmm)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;have a list of all alpha;outputs a log arraybets
(defparameter *hmm-alphabets* nil)

(defmacro def-alphabet (name alphabet)
  `(progn
     (defparameter ,name ,alphabet)
     (push ',name *hmm-alphabets*)))


(def-alphabet *adn-code* '(#\A #\C #\G #\T))
(def-alphabet *arn-code* '(#\A #\C #\G #\U))
(def-alphabet *protein-code* '(#\A #\C #\D #\E #\F #\G #\H #\I #\K #\L #\M #\N #\P #\Q #\R #\S #\T #\V #\W #\Y))

