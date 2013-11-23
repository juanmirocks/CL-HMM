(in-package :cl-hmm)

(defun test-forward&backward_are_symmetric ()
  (semiring :log
    (let* ((N 1)
           (phmm (make-uniform-phmm N 1 1))
           (obs (cbook-encode phmm (list (vector 0) (vector 0))))
           (size_x (length (first obs)))
           (size_y (length (second obs))))
      (multiple-value-bind (o_p alpha) (forward-log phmm obs)
        (multiple-value-bind (beta) (backward-log phmm obs)
          (loop for l to size_x do
               (loop for r to size_y do
                    (when (< 0 (max l r))
                      (loop for i below N
                         with accum = ZERO do
                           (setf accum (log+ accum (log* (alpha[] :log i l r) (beta[] :log i l r))))
                         finally
                           (unless (= accum o_p)
                             (warn "Sum of alpha*beta (~d, ~d) must equal the probability, but ~a != ~a" l r accum o_p))))))
          (values phmm (slot-value phmm 'B) o_p alpha beta))))))

  (test-forward&backward_are_symmetric)
