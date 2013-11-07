(in-package :net.ashrentum.cl-hmm)

(defun hmm-open (filename)
  "Read a hmm file. If only one hmm appears, the output is an hmm object, otherwise is a list of hmm objects"
  (let ((hmms))
    (with-open-file (stream filename)
      (let ((s-expr))
        (while (setf s-expr (read stream nil))
          (push (eval s-expr) hmms))))
    (if (leng=1 hmms) (car hmms) (nreverse hmms))))
