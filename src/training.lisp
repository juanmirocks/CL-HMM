;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Author: Juan Miguel Cejuela
;; Created: Fri Jul 11 23:01:55 2008 (CEST)
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :cl-hmm)

;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Parameter estimation learning algorithms. Currently implemented:
;;   * Baum-Welch (scaled) parameter estimation for HMMs.
;;
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric hmm-estimate-bws
    (hmm obss-c
         &key confidence iterations verbose-estimation
         obss-l starting-noise max-times threshold ri ra rb verbose-bws)
  (:documentation
   "Estimate the model using BaumWelch(-scaled) with the given model parameters and learning observations.
      Repeat the estimation algorithm `iteration`-times starting from a mixture of a random model and the one given (the proportion is balanced with `confidence`).
      hmm: model
      obss-c: list of cbook-encoded observations to train with
      confidence: float in [0, 1], confidence in the given model parameters. 0 to start with a total random model
      iterations: number of iterations
      verbose-estimation:
      obss-l: (optional) list of labeled observations
      starting-noise: float in [0, 1], initial noise for Baum-Welch
      max-times: max number of iterations for Baum-Welch
      threshold: minimum difference with the last loglikelihood to stop the process
      ri: initial probs pseudoconts (vector)
      ra: transition pseudoconts (array)
      rb: emission pseudoconts (array) (If the pseudoconts are not given, these are set to a minimum value not to
      lose any parameter. Set to nil if you do not want this behavior)
      verbose-bws:

      @return (1) estimated model
      @return (2) loglikelihood of the estimated model
      @return (3) list of the observations' loglikelihoods for the estimated model"))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod hmm-estimate-bws
    ((hmm hmm) obss-c
     &key (confidence +es-model-confidence+) (iterations +es-iterations+) (verbose-estimation t)
     obss-l (starting-noise +bw-noise-start+)
     (max-times +bw-max-times+) (threshold +bw-threshold+)
     ri ra rb (verbose-bws nil))
  (let ((best-model)
        (cur-model)
        (best-loglikelihood +most-negative-prob-float+)
        (cur-loglikelihood +most-negative-prob-float+)
        (logs)
        (init-model-noise (prob (- 1 confidence)))
        (bws-iter 0)
        (time0 0))
    (when verbose-estimation
      (format t
              "~2%@@@ Estimate the HMM using Baum-Welch scaled @@@  - i: ~a, c: ~3$, s:~3$, m: ~a, t:~3$~%"
              iterations confidence starting-noise max-times threshold))
    (when (or (< init-model-noise 0) (> init-model-noise 1)) (error "model confidence must be within [0, 1]"))
    (dotimes (i iterations (list best-model best-loglikelihood (nreverse logs)))
      (when verbose-estimation
        (format t "~2%*** iter: # ~a~%" i)
        (format t "---------------------------------------------------------------------------------------------------~%"))
      (setf cur-model (hmm-copy hmm))
      (!hmm-noisify hmm init-model-noise)
      (setq time0 (get-internal-real-time))
      (multiple-value-setq (cur-model cur-loglikelihood bws-iter)
        (funcall (if (eql (type-of hmm) 'phmm) #'baum-welch #'baum-welch-scl) ;trick for now for phmm until -scl is implemented
                 cur-model obss-c :obss-l obss-l
                 :starting-noise starting-noise :max-times max-times :threshold threshold
                 :ri ri :ra ra :rb rb :verbose verbose-bws))
      (when verbose-estimation (format t "~% ->loglikelihood: ~a, no-itr: ~a, time: ~a s~%"
                                       cur-loglikelihood bws-iter (time-elapsed time0)))
      (push cur-loglikelihood logs)
      (when (> cur-loglikelihood best-loglikelihood)
        (setq best-model cur-model
              best-loglikelihood cur-loglikelihood)))))


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Training utilities
;;
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun soften-label-boundaries (prediction wildcard-label &optional (grade 2))
  "Soften the boundaries between the labels regions. Useful when training with labeled sequences of disputable reliability"
    (do ((prediction. (copy-seq prediction))
         (length (length prediction))
         (cur-pos nil next-pos)
         (cur-label (aref prediction 0) (if cur-pos (aref prediction cur-pos)))
         (next-pos 0)
         (sublength 0))
        ((null next-pos) prediction.)
      (setq next-pos (position-if-not #'(lambda (x) (char= x cur-label)) prediction
                                           :start (if cur-pos cur-pos 0))
            sublength (cond
                        ((and cur-pos next-pos) (- next-pos cur-pos))
                        (next-pos next-pos)
                        (cur-pos (- length cur-pos))))
      (loop for g = grade then (1- g)
           until (> sublength (* g 2))
           finally
           (unless (null cur-pos)
             (do ((pos cur-pos (1+ pos))
                  (i 0 (1+ i)))
                 ((= i g))
               (setf (aref prediction. pos) wildcard-label)))
           (unless (null next-pos)
             (do ((pos (- next-pos g) (1+ pos))
                  (i 0 (1+ i)))
                 ((= i g))
               (setf (aref prediction. pos) wildcard-label)))
           (setq cur-pos next-pos))))
