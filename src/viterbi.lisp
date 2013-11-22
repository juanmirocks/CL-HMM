;; Author: Juan Miguel Cejuela
;; Created: Wed Jul  9 19:07:47 2008 (CEST)

(in-package :net.ashrentum.cl-hmm)

;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Viterbi, original and log versions
;;
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric viterbi (hmm obs-c &key labeled ending-with)
  (:documentation
   "Apply the viterbi algorithm to the given observation
	hmm:
	obs-c: observation index-coded (see cbook)
	labeled: NIL (by default), the path uses the internal representation. Otherwise is labeled

	value1: probability of the best path
   	value2: path"))

(defgeneric viterbi-log (hmm obs-c &key labeled ending-with)
  (:documentation
   "Apply the viterbi-log algorithm to the given observation
	hmm:
	obs-c: observation index-coded (see cbook)
	labeled: NIL (by default), the path uses the internal representation. Otherwise is labeled

	value1: log probability of the best path
   	value2: path"))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmacro def-hmm-simple-viterbi-log (&optional log-version)
  `(defmethod ,(if log-version 'viterbi-log 'viterbi)
       ((hmm hmm-simple) obs-c &key labeled ending-with)
     (declare (optimize (speed 3) (safety 0)))
     (restart-case (hmm-incorrect-signal hmm)
       (continue-anyway () nil))
     (hmm-simple-slots (N PE B A iA-to) hmm
       (let* ((PE ,(if log-version `(log-vector PE) 'PE))
              (A ,(if log-version `(log-array A) 'A))
              (B ,(if log-version `(log-array B) 'B))
              (time (length obs-c))
              (phi (make-typed-array (list time N) 'prob-float +most-negative-prob-float+))
              (psi (make-typed-array (list time N) 'cbook-state 0))
              (path (make-typed-array time 'cbook-state 0))
              (last-time (1- time))
              (N-1 (1- N)))
         (declare (cbook-alphabet obs-c) (fixnum time last-time)
                  ((prob-array (* *)) phi A B)
                  ((simple-array cbook-state (* *)) psi)
                  (PE-vec PE))
         ;;Initialisation
         (loop for j fixnum from 0 to N-1 do
              (setf (aref phi 0 j) (,(if log-version '+ '*) (aref PE j) (aref B j (aref obs-c 0)))))
         ;;Induction
         (let ((emi +0-prob+)
               (max +very-negative-prob-float+)
               (new +0-prob+)
               (best-from 0)
               (head 0))
           (declare (prob-float max new emi) (cbook-symbol best-from) (cbook-state head))
           (loop for t0 fixnum = 0 then t+1 for t+1 fixnum from 1 to last-time do
                (loop for j fixnum from 0 to N-1 do
                     (when (> (setq emi (aref B j (aref obs-c t+1))) +very-negative-prob-float+)
                       (setq head (car (aref iA-to j))
                             max (,(if log-version '+ '*) (aref phi t0 head) (aref A head j))
                             best-from head)
                       (dolist-itrans (i (cdr (aref iA-to j)))
                         (setq new (,(if log-version '+ '*) (aref phi t0 i) (aref A i j)))
                         (when (> new max) (setq max new best-from i)))
                       (setf (aref phi t+1 j) (,(if log-version '+ '*) emi max)
                             (aref psi t+1 j) best-from)))))
          ;;Termination
         (multiple-value-bind (prob lpsi) (viterbi-log-end hmm A phi N-1 last-time ending-with)
           (setf (aref path last-time) lpsi)
           ;;Backtracking
           (loop for t+1 = last-time then t0 for t0 from (1- last-time) downto 0 do
                (setf (aref path t0) (aref psi t+1 (aref path t+1))))

           (values (if labeled
                       (look-up path (hmm-state-labels hmm))
                       path)
                   prob))))))

;;;if the hmm is finite it's a step more
(defgeneric viterbi-log-end (hmm-simple A phi N-1 last-time ending-with))

(defmacro def-viterbi-log-end (specializer end-value)
  `(defmethod viterbi-log-end ((hmm ,specializer) A phi N-1 last-time ending-with)
     (declare ((prob-array (* *)) A phi))
     (let ((state-labels (hmm-state-labels hmm)))
       (loop for j fixnum from 0 to N-1
          with max of-type prob-float = +most-negative-prob-float+
          and new of-type prob-float = +0-prob+
          and lpsi fixnum = 0 do
            (when (or (null ending-with) (member (aref state-labels j) ending-with))
              (setf new ,end-value)
              (when (> new max) (setf max new lpsi j)))
          finally (return (values max lpsi))))))

(def-viterbi-log-end hmm-infinite
    (aref phi last-time j))
(def-viterbi-log-end hmm-finite
    (+ (aref phi last-time j) (aref A j N-1))) ;;TODO Flop!! this has to be done as well for the normal version, *!!!

(def-hmm-simple-viterbi-log nil) ;original version
(def-hmm-simple-viterbi-log t) ;log version
