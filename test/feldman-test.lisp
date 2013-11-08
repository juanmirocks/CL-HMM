(in-package :cl-hmm)

(defun protocol-experiment (results-folder &key (states-power 5) (em-num-iterations 6) (max-times 300))
       (multiple-value-bind (in L-size R-size) (read-pair-observations-file "../test/resources/sample_feldman.txt")
         (loop
            for power from 1 to states-power
            for num-states = (expt 2 power)
            do
              (format t "Model with # states: ~i2~%" num-states)
              (let* ((phmm (make-random-phmm num-states L-size R-size))
                     (training-data (mapcar #'(lambda (o) (cbook-encode phmm o)) in))
                     (test-data (mapcar #'(lambda (o) (cbook-encode phmm o)) (read-pair-observations-file "../test/resources/test_feldman.txt"))))
                (multiple-value-bind (best loglikelihood)
                    (hmm-estimate-bws phmm training-data :threshold 0.1 :verbose-bws nil :starting-noise 0 :max-times max-times :iterations em-num-iterations :confidence 0)
                  (list best loglikelihood))))))
