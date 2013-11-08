(in-package :cl-hmm)

(defun protocol-experiment (results-folder &key (states-power 5) (em-num-iterations 6) (max-times 300) (verbose-bws nil))
       (multiple-value-bind (in L-size R-size) (read-pair-observations-file "../test/resources/sample_feldman.txt")
         (loop
            for power from 0 to states-power
            for num-states = (expt 2 power)
            do
              (format t "~3%### Model with num states: ~d~3%" num-states)
              (let* ((phmm (make-random-phmm num-states L-size R-size))
                     (training-data (subseq (mapcar #'(lambda (o) (cbook-encode phmm o)) in) 0 100))
                     (test-data (mapcar #'(lambda (o) (cbook-encode phmm o)) (read-pair-observations-file "../test/resources/test_feldman.txt"))))
                (multiple-value-bind (best loglikelihood)
                    (hmm-estimate-bws phmm training-data :threshold 0.05 :verbose-bws verbose-bws :starting-noise 0 :max-times max-times :iterations em-num-iterations :confidence 0)
                  (with-open-file (stream (format nil "~a/model_~d_~f.txt" results-folder num-states loglikelihood) :direction :output :if-exists :supersede)
                    (loop for obs in test-data
                       for x-encoded = (first obs)
                       for x = (cbook-decode-left best x-encoded)
                       for y = (cbook-decode-right best (hmm-translate best x-encoded))
                       do
                         (loop for i below (1- (length x)) do (format stream "~d_" (elt x i))) (format stream "~d" (elt x (1- (length x))))
                         (format stream " ")
                         (loop for i below (1- (length y)) do (format stream "~d_" (elt y i))) (format stream "~d" (elt y (1- (length y))))
                         (format stream "~%"))))))))
