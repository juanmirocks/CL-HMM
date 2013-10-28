;; Author: Juan Miguel Cejuela
;; Created: Wed Jul  9 19:08:32 2008 (CEST)
;; Last-Updated: 2011-08-11
;;           By: Juan Miguel Cejuela
;;     Update #: 75
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;; Description:
;;
;; Baum-Welch algorithm in the 1) scaled and 2) original versions. By now only
;; for the class hmm inifinite (see hmm-simple). Possible training with labeled
;; sequences. Optional initial noise, with diminishing normalized to the model
;;
;; The code is divided in macros and functions in order to take advantage that
;; both version algorithms are pretty similar, having a common skeleton for
;; both.
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;; Change log:
;; 12-Jul-2008    Ashrentum
;;    Last-Updated: Sat Jul 12 11:09:03 2008 (CEST) #60 (Ashrentum)
;;    Fixed sbcl operation (gc) not compliant with ANSI Common Lisp.
;;    Init Pseudoconts for all parameters no to lost them due to
;;      insufficient training data
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :net.ashrentum.cl-hmm)

(declaim (optimize (speed 3) (safety 0)))
;;(declaim (sb-ext:unmuffle-conditions sb-ext:compiler-note))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(eval-when (:compile-toplevel :load-toplevel :execute)
;;;parameters for the estimation process
  (defconstant +bw-max-times+ 100) ;after baum welch
  (defconstant +bw-threshold+ (prob 0.01))

  ;;noise amplitude, A = S - (i / (log_B C))
  ;;      A, Amplitude
  ;;      S, Starting noise
  ;;      i, iteration
  ;;      C, model Complexity
  ;;      B, noise log Base (to adjust)
  (defconstant +bw-noise-base+ (prob 1.04))
  (defconstant +bw-noise-start+ (prob 0.4))

  (defconstant +es-model-confidence+ (prob 0.7)) ;after estimation
  (defconstant +es-iterations+ 10)

  (defconstant +bw-default-min-pseudocount+ (prob 1d-4)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Baum-Welch, pure&scaled
;;
;; By now only version for HMM infinite are avalaible
;;
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric baum-welch
    (hmm obss-c &key obss-l starting-noise max-times threshold ri ra rb verbose)
  (:documentation "Train the hmm using the pure Baum-Welch algorithm
      hmm: hmm to train
      obss-c: list of cbook-encoded observations to train with
      obss-l: (optional) list of labeled observations
      starting-noise: initial noise to play with (0 to 1)
      max-times: max-times to run the alg.
      threshold: minimum difference change between 2 hmms to accept it and stop
      ri: initial probs pseudoconts (vector)
      ra: transition pseudoconts (array)
      rb: emission pseudoconts (array) (If the pseudoconts are not given, these are set to a minimum value not to
      lose any parameter. Set to nil if you do not want this behavior)
      verbose: prints detailed information about what is happening"))

(defgeneric baum-welch-scl
    (hmm obss-c &key obss-l starting-noise max-times threshold ri ra rb verbose)
  (:documentation "Scaled version of the Baum-Welch algorithm (see baum-welch)"))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;
;;;; The code is divided into common parts to be reused for both the original
;;;; and the scaled versions of the Baum-Welch algorithm
;;;;

;;;;variable definition (not very functional indeed)
(defmacro hmm-simple-vars (scaled &body body)
  `(hmm-simple-slots (PE N M no-groups state-groups A B iA-from) hmm
     (let ((nPE (make-typed-array N 'prob-float +0-prob+)) ;named after new PE
           (PErow +0-prob+) ;(see 3.18 Durbin BSA)
           (nA (make-typed-array (list N N) 'prob-float +0-prob+))
           (Arow (make-typed-array N 'prob-float +0-prob+))
           (nB (make-typed-array (list N M) 'prob-float +0-prob+))
           (Brow (make-typed-array N 'prob-float +0-prob+))
           (last-loglikelihood +most-negative-prob-float+)
           (cur-loglikelihood +most-negative-prob-float+)
           (x^k_size 0)
           ;;pseudoconts. If not given, set them an uniform value not to lose any parameter due to insufficient training
           (ri (cond
                 (ri ri)
                 ((and (not ri) rip) nil)
                 (t (make-typed-array N 'prob-float +bw-default-min-pseudocount+))))
           (ra (cond
                 (ra ra)
                 ((and (not ra) rap) nil)
                 (t (make-typed-array (list N N) 'prob-float +bw-default-min-pseudocount+))))
           (rb (cond
                 (rb rb)
                 ((and (not rb) rbp) nil)
                 (t (make-typed-array (list N M) 'prob-float +bw-default-min-pseudocount+))))
           (time0 0) ; to measure the algorithm's running time
           (noise-amp (prob starting-noise)) ;noise amplitude
           (noise +0-prob+)
           (noise-decrease (prob (/ (log (hmm-complexity hmm) +bw-noise-base+))))
           ,@(if scaled `((scale (make-typed-array '(0) 'prob-float +0-prob+)) (P{x^k} +0-prob+)) `((1/P{x^k} +0-prob+))))
       (declare ((prob-array (*)) nPE Arow Brow)
                ((prob-array (* *)) nA nB)
                (prob-float PErow last-loglikelihood cur-loglikelihood noise-amp noise-decrease)
                (fixnum x^k_size time0))
       (declare ,@(if scaled `(((prob-array (*)) scale) (prob-float P{x^k})) `((prob-float 1/P{x^k}))))
       ,@body)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(eval-when (:compile-toplevel :load-toplevel :execute)
;;; Add the pseudocounts if not nil
  (defun hmm-simple-pseudoconts ()
    `((when ri (dotimes (i N)
                 (setf (aref (the PE-vec nPE) i) (aref ri i))
                 (incf PErow (aref ri i))))
      (when ra (dotimes (i N)
                 (dotimes (j N)
                   (setf (aref nA i j) (aref ra i j))
                   (incf (aref Arow i) (aref ra i j)))))
      (when rb (dotimes (i N)
                 (dotimes (b M)
                   (setf (aref nB i b) (aref rb i b))
                   (incf (aref Brow i) (aref rb i b)))))))

;;; Group the emissions probabilities
  (defun group-emissions (N M no-groups state-groups nB Brow)
    (declare (cbook-state N M no-groups) (B-array nB) ((prob-array (*)) Brow))
    (let* ((nB-grouped (make-typed-array (list no-groups M) 'prob-float +0-prob+))
           (Brow-grouped (make-typed-array no-groups 'prob-float +0-prob+)))
      (dotimes (i N)
        (let ((group-i (aref state-groups i)))
          (declare (cbook-state group-i))
          (dotimes (j M)
            (incf (aref nB-grouped group-i j) (aref nB i j)))
          (incf (aref Brow-grouped group-i) (aref Brow i))))
      (values nB-grouped Brow-grouped)))

;;; Update parameters
  (defun hmm-simple-update ()
    `(multiple-value-bind (nB-grouped Brow-grouped)
         (group-emissions N M no-groups state-groups nB Brow)
       (do* ((i 0 (1+ i))
             (pArow (aref Arow i) (aref Arow i))
             (gi (aref state-groups i) (aref state-groups i))
             (pBrow-grouped (aref Brow-grouped gi) (aref Brow-grouped gi)))
            ((= i N) nil)
         (declare (cbook-state i) (prob-float pArow pBrow-grouped))
         ;;init
         (if (= +0-prob+ PErow) (error "No info for init probabilities!")
             (setf (aref PE i) (/ (aref nPE i) PErow)
                   (aref nPE i) +0-prob+)) ;reset
         ;;transitions
         (if (= +0-prob+ pArow) (error "No info for trans probability in state number: ~a" i)
             (dolist-itrans (j (aref iA-from i))
               (setf (aref A i j) (/ (aref nA i j) pArow)
                     (aref nA i j) +0-prob+))) ;reset
         ;;emissions GROUPED
         (if (= +0-prob+ pBrow-grouped) (error "No info for emis probability in state number: ~a" i)
             (dotimes (s M)
               (setf (aref B i s) (/ (aref nB-grouped gi s) pBrow-grouped)
                     (aref nB i s) +0-prob+))) ;reset
         ;;reset
         (setf (aref Arow i) +0-prob+)
         (setf (aref Brow i) +0-prob+))
       (setf PErow +0-prob+))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Define a Baum-Welch method
(defmacro def-baum-welch (name specializer-hmm scaled &body algorithm)
  `(defmethod ,name
       ((hmm ,specializer-hmm)
        obss-c
        &key
          obss-l
          (starting-noise +bw-noise-start+)
          (max-times +bw-max-times+) (threshold +bw-threshold+)
          (ri nil rip) (ra nil rap) (rb nil rbp)
          (verbose nil))

     (declare (optimize (speed 3) (safety 0)) (fixnum max-times) (float threshold) (list obss-c))
     (restart-case (hmm-incorrect-signal hmm)
       (continue-anyway () nil))
     (setf hmm (hmm-copy hmm)) ;don't overwrite the given hmm
     (hmm-simple-vars ,scaled

       ;;Recurrence. Termination when max-times or likelihood is less than the threshold
       (loop for z fixnum from 1 to max-times do
            (setq time0 (get-internal-real-time))
            (setq last-loglikelihood cur-loglikelihood
                  cur-loglikelihood +0-prob+)

            ,@(hmm-simple-pseudoconts) ;constant pseudoconts

            (do ((obss obss-c (cdr obss))
                 (x^k (make-typed-array 0 'cbook-symbol 0))
                 (obssl obss-l (cdr obssl))
                 (x^k-labels (make-typed-array 0 'state-label +label-wildcard+)))
                ((null obss) nil)
              (declare (cbook-alphabet x^k))
              (setf x^k (car obss)
                    x^k_size (length x^k)
                    x^k-labels (car obssl))

              ;; ----------------------------------------------------------------------
              ,@algorithm
              ;; ----------------------------------------------------------------------

              #+sbcl (sb-ext:gc :gen 1 :full t)) ;free memory

            ,(hmm-simple-update) ;update parameters

          ;; Apply noise
            (setf noise (* (random +1-prob+) noise-amp))
            (!hmm-noisify hmm noise)
            (decf noise-amp noise-decrease)
            (when (< noise-amp 0) (setf noise-amp +0-prob+))

            (when verbose
              (format t "~a:~5T~a~28T noise: ~3$  (~3$ s)"
                      z cur-loglikelihood noise (time-elapsed time0))
              (when (< cur-loglikelihood last-loglikelihood)
                (format t "   worse! (~a)" (- cur-loglikelihood last-loglikelihood)))
              (fresh-line))

          until (or
                 (and (< (abs (- cur-loglikelihood last-loglikelihood)) threshold)
                      (> cur-loglikelihood last-loglikelihood)
                      (zerop noise-amp))
                 (zerop cur-loglikelihood))

          finally
            (multiple-value-bind (correct details) (hmm-correctp hmm)
              (unless correct
                (warn "Some states are dumb now. Output of hmm-correct-p:~2%~a~%" details)))
            (return (values hmm cur-loglikelihood z))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; the actual core of the algorithm (2 versions, pure & scaled)
(defmacro hmm-infinite-algorithm-core (&key (scaled nil))
  ;;for all the states
  `(loop for i of-type cbook-state from 0 below N do
      ;;init distribution
        (let ((pPErow (* (aref alphas 0 i) (aref betas 0 i))))
          ,(if scaled
               `(if (/= +0-prob+ (aref scale 0))
                    (setq pPErow (/ pPErow (aref scale 0)))
                    (error "Error, value for scale 0 is 0. Check the model or the labels"))
               `(*= pPErow 1/P{x^k}))
          (incf (aref nPE i) pPErow)
          (incf PErow pPErow))
      ;;transitions
        (if (leng=1 (aref iA-from i))
            (progn
              (incf (aref nA i (car (aref iA-from i))) +1-prob+) ;;the probability is fixed, it's always 1
              (incf (aref Arow i) +1-prob+))
            (dolist-itrans (j (aref iA-from i))
              (loop for t0 = 0 then t+1 for t+1 from 1 below x^k_size
                 sum (* (aref alphas t0 i)
                        (aref B j (aref x^k t+1))
                        (aref betas t+1 j)) into pArow of-type prob-float
                 finally
                   ,(if scaled
                        `(*= pArow (aref A i j))
                        `(*= pArow (* (aref A i j) 1/P{x^k})))
                   (incf (aref nA i j) pArow)
                   (incf (aref Arow i) pArow))))
      ;;emissions
        (loop for s of-type cbook-state from 0 below M for emis of-type prob-float = (aref B i s) then (aref B i s) do
             (unless (zerop emis)
               (loop for t0 from 0 below x^k_size with pBrow of-type prob-float = +0-prob+ do
                    (when (= s (aref x^k t0))
                      ,(if scaled
                           `(if (/= +0-prob+ (aref scale t0))
                                (incf pBrow (/ (* (aref alphas t0 i) (aref betas t0 i)) (aref scale t0)))
                                (error "Error, value for scale ~d is 0. Check the model or the labels" t0))
                           `(incf pBrow (* (aref alphas t0 i) (aref betas t0 i)))))
                  finally
                    ,@(if scaled
                          `((incf (aref nB i s) pBrow))
                          `((*= pBrow 1/P{x^k})
                            (incf (aref nB i s) pBrow)))
                    (incf (aref bRow i) pBrow))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;baumm-welch pure :: hmm-infinite
(def-baum-welch baum-welch hmm-infinite nil
  (let ((alphas (make-typed-array '(0 0) 'prob-float +0-prob+))
        (betas (make-typed-array '(0 0) 'prob-float +0-prob+)))
    (declare ((prob-array (* *)) alphas betas)) ;we declare here the alphas and betas for memory issues
    (multiple-value-setq (1/P{x^k} alphas) (forward hmm x^k))
    (setf betas (backward hmm x^k))
    (unless (zerop 1/P{x^k})
      (incf cur-loglikelihood (the prob-float (log 1/P{x^k})))
      (setf 1/P{x^k} (/ 1/P{x^k})) ;now yes, it's the inverse
      (hmm-infinite-algorithm-core :scaled nil))))

;;baum-welch scaled :: hmm-infinite
(def-baum-welch baum-welch-scl hmm-infinite t
  (let ((alphas (make-typed-array '(0 0) 'prob-float +0-prob+))
        (betas (make-typed-array '(0 0) 'prob-float +0-prob+)))
    (declare ((prob-array (* *)) alphas betas))
    (multiple-value-setq (P{x^k} alphas scale) (forward-scl hmm x^k x^k-labels))
    (setf betas (backward-scl hmm x^k scale x^k-labels))
    (incf cur-loglikelihood P{x^k})
    (hmm-infinite-algorithm-core :scaled t)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun pseudocounts (var varp size)
  (cond
    (var var)
    (varp nil)
    (t (make-typed-array size 'prob-float +bw-default-min-pseudocount+))))

(defmethod baum-welch
    ((hmm phmm)
     obss-c
     &key
       obss-l
       (starting-noise +bw-noise-start+)
       (max-times +bw-max-times+) (threshold +bw-threshold+)
       (ri nil rip) (ra nil rap) (rb nil rbp)
       (verbose nil))

  (setf hmm (hmm-copy hmm)) ;don't overwrite the given hmm

  (phmm-slots (N L-size R-size) hmm
    (let ((ri (pseudocounts ri rip N))
          (ra (pseudocounts ra rap `(,N ,N)))
          (rb (pseudocounts rb rbp `(,N ,L-size ,R-size))))

    (loop for iteration from 1
       for time-itr-start = (get-internal-real-time)
       for last-loglikelihood = +most-negative-prob-float+ then cur-loglikelihood
       for cur-loglikelihood = +most-negative-prob-float+
       for noise-amplitude = (max 0 (- starting-noise (/ iteration (log (hmm-complexity hmm) +bw-noise-base+))))
       for noise = (* (random +1-prob+) noise-amplitude)
       do
         (loop for o in obss-c
            for x = (first o)
            for y = (second o)
            for size_x fixnum = (length x)
            for size_y fixnum = (length y)
            for (o_likelihood alpha) = (multiple-value-list (forward hmm o))
            for beta = (backward hmm o)
            do

              (print o))

         (when verbose
           (format t "~a:~5T~a~28T noise: ~3$  (~3$ s)" iteration cur-loglikelihood noise (time-elapsed time-itr-start))
           (when (< cur-loglikelihood last-loglikelihood)
             (format t "   worse! (~a)" (- cur-loglikelihood last-loglikelihood)))
           (fresh-line))

       until (or (= iteration max-times)
                 (and
                  (< (abs (- cur-loglikelihood last-loglikelihood)) threshold)
                  (> cur-loglikelihood last-loglikelihood))
                 (zerop cur-loglikelihood)) ;perfect model to the training data

       finally
         (multiple-value-bind (correct details) (hmm-correctp hmm)
           (unless correct (warn "The output model is incorrect. Output of hmm-correct-p:~2%~a~%" details)))
         (return (values hmm cur-loglikelihood iteration))))))
