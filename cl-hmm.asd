;; Author: Juan Miguel Cejuela
;; Created: Wed Jul  9 19:20:39 2008 (CEST)
;; Last-Updated: 2011-08-11
;;           By: Juan Miguel Cejuela
;;     Update #: 37

(defpackage :net.ashrentum.cl-hmm-system (:use :asdf :cl))
(in-package :net.ashrentum.cl-hmm-system)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defsystem cl-hmm
    :name "cl-hmm"
    :author "Juan Miguel Cejuela"
    :version "0.2"
    :maintainer "Juan Miguel Cejuela"
    :licence "LGPL 3"
    :description "Simple HMM library for Common Lisp"
    :components
    ((:module
      "src"
      :components
      ((:file "packages")
       (:file "cl-hmm"
        :depends-on ("packages"))
       (:file "utilities"
        :depends-on ("packages" "cl-hmm"))

       (:file "hmm-simple"
        :depends-on ("packages" "cl-hmm" "utilities"))

       (:file "viterbi"
        :depends-on ("packages" "cl-hmm" "utilities" "hmm-simple"))
       (:file "for&back-ward"
        :depends-on ("packages" "cl-hmm" "utilities" "hmm-simple"))
       (:file "baum-welch"
        :depends-on ("packages" "cl-hmm" "utilities" "hmm-simple" "for&back-ward"))
       (:file "training"
        :depends-on ("packages" "cl-hmm" "utilities" "hmm-simple" "for&back-ward" "baum-welch"))

       (:file "alphabets"
        :depends-on ("packages"))
       (:file "hmm-files"
        :depends-on ("packages" "cl-hmm" "hmm-simple" "alphabets")))))

    :depends-on
    (:jmc.cl.utils))

