(in-package :cl-hmm)

(defun levenshtein (str1 str2)
  "Calculates the Levenshtein distance between str1 and str2, returns an editing distance (int)."
  (declare (optimize (speed 3) (safety 0)))
  (let ((n (length str1))
        (m (length str2)))
    ;; Check trivial cases
    (cond ((= 0 n) (return-from levenshtein m))
          ((= 0 m) (return-from levenshtein n)))
    (let ((col (make-svector (1+ m)))
          (prev-col (make-svector (1+ m))))
      ;; We need to store only two columns---the current one that
      ;; is being built and the previous one
      (dotimes (i (1+ m))
        (setf (svref prev-col i) i))
      ;; Loop across all chars of each string
      (dotimes (i n)
        (setf (svref col 0) (1+ i))
        (dotimes (j m)
          (setf (svref col (1+ j))
                (min (1+ (svref col j))
                     (1+ (svref prev-col (1+ j)))
                     (+ (svref prev-col j)
                        (if (eq (aref str1 i) (aref str2 j)) 0 1)))))
        (rotatef col prev-col))
      (the fixnum (svref prev-col m)))))

;; A few simple test-cases
(assert (zerop (levenshtein "kitten" "kitten")))
(assert (= (levenshtein "kitten" "") 6))
(assert (= (levenshtein "kitten" "sitting") 3))
(assert (= (levenshtein "" "") 0))
(assert (= (levenshtein "" "a") 1))
(assert (= (levenshtein "a" "") 1))

(defun wer (reference recognized)
  (declare (optimize (speed 3) (safety 0)) (inline wer))
  (min 1.0 (coerce (/ (levenshtein reference recognized) (length reference)) 'float)))

(defun wer-phmm (phmm test-data num-translations)
  (loop for obs in test-data
     with wer = 0.0
     with test-size = (length test-data)
     do
       (loop repeat num-translations
          with o_wer = 0.0
          for x = (first obs)
          for y = (second obs)
          for y-translation = (hmm-translate phmm x)
          do
            (incf o_wer (wer y y-translation))
          finally
            (incf wer (/ o_wer num-translations)))
     finally (return  (/ wer test-size))))

(defun protocol-experiment (results-folder &key (training-size 10000) (num-translations 1000) (states-power 5) (em-num-iterations 6) (max-times 300) (verbose-bws nil))
       (multiple-value-bind (in L-size R-size) (read-pair-observations-file (pwd "test/resources/sample_feldman.txt"))
         (loop
            for power from 0 to states-power
            for num-states = (expt 2 power)
            do
              (format t "~3%### Model with num states: ~d~3%" num-states)
              (let* ((phmm (make-random-phmm num-states L-size R-size))
                     (training-data (subseq (mapcar #'(lambda (o) (cbook-encode phmm o)) in) 0 training-size))
                     (test-data (mapcar #'(lambda (o) (cbook-encode phmm o)) (read-pair-observations-file (pwd "test/resources/test_feldman.txt")))))
                (multiple-value-bind (best loglikelihood)
                    (hmm-estimate-bws phmm training-data :threshold 0.05 :verbose-bws verbose-bws :starting-noise 0 :max-times max-times :iterations em-num-iterations :confidence 0)
                  (hmm-save best (format nil "~a/model_~d_~d_~f.lisp" results-folder num-states training-size loglikelihood) :complete)
                  (with-open-file (stream (format nil "~a/wer_~d_~d_~f.txt" results-folder num-states training-size loglikelihood) :direction :output :if-exists :supersede)
                    (format stream "~f~%" (wer-phmm best test-data num-translations))))))))
