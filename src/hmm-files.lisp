;; Author: Juan Miguel Cejuela
;; Maintainer:
;; Created: Wed Jul  9 19:09:34 2008 (CEST)
;; Last-Updated: 2011-08-11
;;           By: Juan Miguel Cejuela
;;     Update #: 4

(in-package :net.ashrentum.cl-hmm)
(declaim (optimize (speed 3) (safety 0)))


(defmethod hmm-save (filename (hmm hmm-simple) &optional (model-spec :relevant))
  (labels ((model-spec-relevant (hmm states)
             (hmm-simple-slots (V S PE A B) hmm
               (let ((trans (multiple-value-bind (a b) (trans-stats A S nil t) a b))
                     (emis (emis-stats B V S t))
                     (model-spec))
                 (loop
                    for i = 0 then (1+ i)
                    for s in states
                    for a in trans
                    for e in emis
                    for ia = nil then nil
                    for ie = nil then nil do
                      (loop for a. in (cddr a) do
                           (push a. ia))
                      (unless (null (cdddr e))
                        (loop for e. in (cddr e) do
                             (when (numberp e.)
                               (push e. ie))))
                      (if ie
                          (push (list (first s) (second s) (aref PE i) (nreverse ia) (nreverse ie)) model-spec)
                          (push (list (first s) (second s) (aref PE i) (nreverse ia)) model-spec)))
                 (nreverse model-spec)))))

    (unless (typep hmm 'hmm)
      (error "The object is not a valid hmm"))
    (with-open-file (stream filename :direction :output :if-does-not-exist :create :if-exists :rename)
      (prin1
       (ecase (type-of hmm)
         ((hmm-infinite hmm-finite hmm-simple)
          (hmm-simple-slots (S PE A B) hmm
            (let ((states))
              (dotimes (i (hmm-no-states hmm) (setq states (nreverse states)))
                (push (list (state-name (aref S i)) (state-label (aref S i))) states))

              `(make-hmm-simple
                ,(hmm-no-states hmm)
                ,(hmm-no-emissions hmm)
                ',(sequence->list (hmm-alphabet hmm))
                ,(ecase model-spec
                        (:complete
                         `'(,states
                            ,(sequence->list PE)
                            ,(matrix->list A)
                            ,(matrix->list B)))
                        (:relevant `',(model-spec-relevant hmm states)))

                :name ,(hmm-name hmm)
                :alphabet-type ',(array-element-type (hmm-alphabet hmm))
                :model-spec ,model-spec)))))
       stream)) T))


(defun hmm-open (filename)
  "Read a hmm file. If only one hmm appears, the output is an hmm object, otherwise is a list of hmm objects"
  (let ((hmms))
    (with-open-file (stream filename)
      (let ((s-expr))
        (while (setf s-expr (read stream nil))
          (push (eval s-expr) hmms))))
    (if (leng=1 hmms) (car hmms) (nreverse hmms))))
