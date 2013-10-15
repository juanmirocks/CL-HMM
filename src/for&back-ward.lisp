;; Author: Juan Miguel Cejuela
;; Created: Wed Jul  9 19:06:53 2008 (CEST)
;; Last-Updated: 2011-08-11
;;           By: Juan Miguel Cejuela
;;     Update #: 32

(in-package :net.ashrentum.cl-hmm)
;;(declaim (sb-ext:unmuffle-conditions sb-ext:compiler-note))


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Forward and Backward (original)
;;
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric forward (hmm obs-c)
  (:documentation
   "Forward algorithm to the observation
	hmm
	obs-c: observation index-coded (see cbook)
	value1: probability
	value2: alphas matrix"))

;;simple wrapper by now. Possible to don't calculate all alpha values to get this value faster
(defgeneric likelihood (hmm obs-c)
  (:documentation
   "Likelihood to the observation
	hmm:
	obs-c: observation index-coded
	value: likehihood"))

(defgeneric backward (hmm obs-c)
  (:documentation
   "Backward algorithm to the observation
	hmm: model
	obs-c: observation index-coded (see cbook)
	value: betas matrix"))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod forward ((hmm hmm-simple) obs-c)
  (declare (optimize (speed 3) (safety 0)) (cbook-alphabet obs-c))
  (hmm-simple-slots (N PE A B iA-to) hmm
    (let* ((time (length obs-c))
           (alpha (make-typed-array (list time N) 'prob-float +0-prob+))
           (N-1 (1- N))
           (last-time (1- time)))
      (declare (fixnum N time N-1 last-time))
      ;;Initialisation
      (loop for i from 0 to N-1 do
           (setf (aref alpha 0 i) (prob (* (aref PE i) (aref B i (aref obs-c 0))))))
      ;;Induction
      (let ((emi +0-prob+)
            (sum +0-prob+))
        (declare (prob-float emi sum))
        (loop for t0 fixnum = 0 then t+1 for t+1 fixnum from 1 to last-time do
             (loop for j fixnum from 0 to N-1 do
                  (when (> (setq emi (aref B j (aref obs-c t+1))) +0-prob+)
                    (setq sum +0-prob+)
                    (dolist-itrans (i (aref iA-to j))
                      (setq sum (+ sum (* (aref alpha t0 i) (aref A i j)))))
                    (setf (aref alpha t+1 j) (* emi sum))))))
      ;;Termination
      (let ((prob (fin-forward hmm A alpha N-1 last-time)))
        (values prob alpha)))))

;;submethods
(defgeneric fin-forward (hmm-simple A alpha N-1 last-time))
(defmethod fin-forward ((hmm hmm-infinite) A alpha N-1 last-time)
  (declare ((prob-array (* *)) alpha A) (fixnum N-1 last-time))
  (loop for i fixnum from 0 to N-1
     sum (aref alpha last-time i) into sum of-type prob-float
     finally (return sum)))
(defmethod fin-forward ((hmm hmm-finite) A alpha N-1 last-time)
  (declare ((prob-array (* *)) alpha A) (fixnum N-1 last-time))
  (loop for i fixnum from 0 to (1- N-1) ;must finish
     sum (* (aref A i N-1) (aref alpha last-time i)) into sum of-type prob-float
     finally (return sum)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod likelihood ((hmm hmm-simple) obs-c)
  (declare (optimize (speed 3) (safety 0)) (cbook-alphabet obs-c))
  (values (forward hmm obs-c)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod backward ((hmm hmm-simple) obs-c)
  (declare (optimize (speed 3) (safety 0)) (cbook-alphabet obs-c))
  (hmm-simple-slots (N A B iA-from) hmm
    (let* ((time (length obs-c))
           (beta (make-typed-array (list time N) 'prob-float +0-prob+))
           (N-1 (1- N))
           (last-time (1- time)))
      (declare (fixnum time N-1 last-time) ((prob-array (* *)) beta))
      ;;Initialization
      (setf beta (init-backward hmm A beta N-1 last-time))
      ;;Induction
      (let ((sum +0-prob+))
        (declare (prob-float sum))
        (loop for t+1 fixnum = last-time then t0 for t0 fixnum from (1- last-time) downto 0 do
             (loop for i fixnum from 0 to N-1 do
                  (setq sum +0-prob+)
                  (dolist-itrans (j (aref iA-from i))
                    (setq sum (+ sum (* (aref A i j)
                                        (aref b j (aref obs-c t+1))
                                        (aref beta t+1 j)))))
                (loop for j fixnum from 0 to N-1
                   sum (* (aref A i j) (* (aref B j (aref obs-c t+1)) (aref beta t+1 j))) into sum of-type prob-float
                   finally (setf (aref beta t0 i) sum))))
        beta))))

;;submethods
(defgeneric init-backward (hmm-simple A beta N-1 last-time))
(defmethod init-backward ((hmm hmm-infinite) A beta N-1 last-time)
  (declare ((prob-array (* *)) beta A) (fixnum N-1 last-time))
  (loop for j fixnum from 0 to N-1 with iprob = +1-prob+ do
       (setf (aref beta last-time j) iprob))
  beta)
(defmethod init-backward ((hmm hmm-finite) A beta N-1 last-time)
  (declare ((prob-array (* *)) beta A) (fixnum N-1 last-time))
  (loop for j fixnum from 0 to N-1 do ;must finish
       (setf (aref beta last-time j) (aref A j N-1)))
  beta)


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Forward and Backward (scaled)
;;
;; for this moment are only written the infinite versions
;;
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric forward-scl (hmm obs-c &optional obs-l)
  (:documentation
   "Apply the forward scaled algorithm to the given observation
	hmm: model
	obs-c: observation index-coded (see cbook)
	obs-l: labeled observation
	value1: log probability
	value2: alpha matrix
	value3: scale"))

;;simple wrapper by now. Possible to don't calculate all alpha values to get this value faster
(defgeneric loglikelihood (hmm obs-c)
  (:documentation
   "Loglikelihood of the observation for the given model
	hmm: model
	obs-c: observation index-coded
	obs-l: labeled observation
	value: loglikehihood"))

(defgeneric backward-scl (hmm obs-c scale &optional obs-l)
  (:documentation
   "Apply the forward scaled algorithm to the given observation
	hmm: model
	obs-c: observation index-coded (see cbook)
	scale: outputed by the forward-sc
	obs-l: labeled observation
	value: beta matrix"))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod forward-scl ((hmm hmm-simple) obs-c &optional obs-l)
  (declare (optimize (speed 3) (safety 0)) (cbook-alphabet obs-c))
  (restart-case (hmm-incorrect-signal hmm)
    (continue-anyway () nil))
  (hmm-simple-slots (N PE A iA-to B) hmm
    (let* ((time (length obs-c))
           (alpha-scl (make-typed-array (list time N) 'prob-float +0-prob+))
           (scale (make-typed-array time 'prob-float +0-prob+)) ;scaling factor vector
           (last-time (1- time))
           (N-1 (1- N))
           (state-labels (hmm-state-labels hmm)))
      (declare (cbook-alphabet obs-c) (fixnum time last-time) ((simple-array state-label (*)) state-labels))
      ;;Initialisation
      ;;(setf (aref scale last-time) +1-prob+)
      (loop for j from 0 to (the fixnum N-1) do
           (when (or (null obs-l) (compatible-label (aref state-labels j) (aref (the string obs-l) 0)))
             (setf (aref alpha-scl 0 j) (* (aref PE j) (aref B j (aref obs-c 0)))))
           (incf (aref scale 0) (aref alpha-scl 0 j))) ;factor scale
      (if (= +0-prob+ (aref scale 0))
          (return-from forward-scl (values +very-negative-prob-float+ alpha-scl scale));;(setf (aref scale 0) +1-prob+)
          (setf (aref scale 0) (/ (aref scale 0))))
      (loop for j from 0 to (the fixnum N-1) with factor = (aref scale 0) do
           (*= (aref alpha-scl 0 j) factor)) ;normalized
      ;;Induction
      (let ((emi +0-prob+)
            (sum +0-prob+))
        (declare (prob-float emi sum))
        (loop for t0 = 0 then t+1 for t+1 from 1 to last-time do
             (loop for j fixnum from 0 to N-1 do
                  (when (and (> (setq emi (aref b j (aref obs-c t+1))) +very-negative-prob-float+)
                             (or (null obs-l) (compatible-label (aref state-labels j) (aref (the string obs-l) t+1))))
                    (setq sum +0-prob+)
                    (dolist-itrans (i (aref iA-to j))
                      (incf sum (* (aref alpha-scl t0 i) (aref A i j))))
                    (setf (aref alpha-scl t+1 j) (* sum emi))
                    (incf (aref scale t+1) (aref alpha-scl t+1 j))))
             (when (= +0-prob+ (aref scale t+1))
                 (warn "The value for scale in ~d is 0. Check the model or the labels" t+1)
                 (return-from forward-scl (values +most-negative-prob-float+ alpha-scl scale)))
             (setf (aref scale t+1) (/ (aref scale t+1)))
             (loop for j from 0 to (the fixnum N-1) with factor = (aref scale t+1) do
                  (*= (aref alpha-scl t+1 j) factor))))
      ;;Termination
      (let ((prob (loop for t0 fixnum from 0 to last-time
                     sum (the prob-float (log (aref scale t0))) into sum of-type prob-float
                     finally (return (- sum)))))
        (values prob alpha-scl scale)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod loglikelihood ((hmm hmm-simple) obs-c)
  (declare (optimize (speed 3) (safety 0)) (cbook-alphabet obs-c))
  (values (forward-scl hmm obs-c)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod backward-scl ((hmm hmm-simple) obs-c scale &optional obs-l)
  (declare (optimize (speed 3) (safety 0)) (cbook-alphabet obs-c) ((prob-array (*)) scale))
  (restart-case (hmm-incorrect-signal hmm)
    (continue-anyway () nil))
  (hmm-simple-slots (N A iA-from B) hmm
    (let* ((time (length obs-c))
           (last-time (1- time))
           (beta-scl (make-typed-array (list time N) 'prob-float +0-prob+))
           (N-1 (1- N))
           (state-labels (hmm-state-labels hmm)))
      (declare (fixnum time last-time) (cbook-state N-1) (itrans iA-from) ((simple-array state-label (*)) state-labels))
      ;;Initialization
      (loop for i from 0 to (the fixnum N-1) do
           (when (or (null obs-l) (compatible-label (aref state-labels i) (aref (the simple-string obs-l) last-time)))
             (setf (aref beta-scl last-time i) (aref scale last-time))))
      ;;Induction
      (let ((sum +0-prob+))
        (declare (prob-float sum))
        (loop for t+1 = last-time then t0 for t0 from (1- last-time) downto 0 do
             (loop for i fixnum from 0 to N-1 do
                  (setq sum +0-prob+)
                  (when (or (null obs-l) (compatible-label (aref state-labels i) (aref (the simple-string obs-l) t0)))
                    (dolist-itrans (j (aref iA-from i))
                      (setq sum (+ sum (* (aref A i j)
                                          (aref b j (aref obs-c t+1))
                                          (aref beta-scl t+1 j)))))
                    (setf (aref beta-scl t0 i) (* sum (aref scale t0)))))
           finally (return beta-scl))))))


;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Posterior Probabilities
;;
;; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric posterior-probs-labels (hmm obs-c &optional ending-with obs-l)
  (:documentation
   "Give the posterior probablities grouped by the labels"))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;use ending-with to forward and backward?
(defmethod posterior-probs-labels ((hmm hmm-simple) obs-c &optional ending-with obs-l)
  (multiple-value-bind (a b) (hmm-correctp hmm)
    (unless a (error 'hmm-incorrect :text b)))
  (let* ((N (hmm-no-states hmm))
         (time (length obs-c))
         (state-labels (hmm-state-labels hmm))
         (state-labels-indexed (make-typed-array N 'cbook-state 0))
         (group-labels (delete-duplicates state-labels :from-end t))
         (no-labels (length group-labels))
         (labels-hash (make-hash-table :size no-labels :test 'equalp))
         (probs-label-grouped (make-typed-array (list time no-labels) 'prob-float +0-prob+))
         (path (make-typed-array time 'state-label #\.)))
    (declare (cbook-alphabet obs-c) (fixnum time) ((simple-array state-label (*)) state-labels))
    (dotimes (i no-labels) ;get the label groups
      (setf (gethash (aref group-labels i) labels-hash) i))
    (dotimes (i N) ;set the label group index for every state
      (setf (aref state-labels-indexed i) (gethash (aref state-labels i) labels-hash)))
    ;;put together the posterior probablities for the label groups
    (multiple-value-bind (loglikelihood alpha-scl scale) (forward-scl hmm obs-c obs-l)
      (multiple-value-bind (beta-scl) (backward-scl hmm obs-c scale obs-l)
        (dotimes (t0 time)
          (dotimes (i N)
             (incf (aref probs-label-grouped t0 (aref state-labels-indexed i))
                   (/ (* (aref alpha-scl t0 i) (aref beta-scl t0 i)) (aref scale t0))))
          (let ((max +very-negative-prob-float+)
                (from #\.))
            (dotimes (l no-labels)
              (when (> (aref probs-label-grouped t0 l) max)
                (setf max (aref probs-label-grouped t0 l)
                      from l)))
            (setf (aref path t0) (aref group-labels from)))))
     (values group-labels probs-label-grouped path loglikelihood))))








