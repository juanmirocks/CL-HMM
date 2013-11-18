(in-package :cl-hmm)

(defun forward_performance (&key (states-power 5) (size-data 10000))
  (multiple-value-bind (in L-size R-size) (read-pair-observations-file (pwd "test/resources/sample_feldman.txt"))
    (time
     (loop
        for power from 0 to states-power
        for num-states = (expt 2 power)
        do
          (format t "### Model with num states: ~d~%" num-states)
          (let* ((phmm (make-random-phmm num-states L-size R-size))
                 (training-data (subseq (mapcar #'(lambda (o) (cbook-encode phmm o)) in) 0 size-data)))
            (loop for o in training-data do
                 (forward phmm o)))))))
