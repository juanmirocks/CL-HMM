;; Author: Juan Miguel Cejuela
;; Created: Wed Jul  9 19:13:54 2008 (CEST)
;; Last-Updated: 2011-08-11
;;           By: Juan Miguel Cejuela
;;     Update #: 15

(in-package :cl-hmm)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;fastest search in a prob-float array, (see jmc.cl.utils)
(def-lin-search array-search prob-float (simple-array prob-float) >=)

(defun select-random (accum-array &key (indices-1 nil) (fixed-max nil) (fixed-pick nil))
  (declare (optimize (speed 3) (safety 0)))
  (let ((dims (array-dimensions accum-array)))
    (unless (= (length dims) (1+ (length indices-1))) (error "The size of indices-1 must be exactly array's dimensions - 1"))
    (labels ((fstart (indices-1 dims-1 accum)
               (if (null indices-1)
                   accum
                   (fstart (cdr indices-1) (cdr dims-1) (+ accum (apply #'* (car indices-1) dims-1))))))
      (let* ((start (fstart indices-1 (cdr dims) 0))
             (end (+ start (car (last dims))))
             (max (if fixed-max fixed-max (row-major-aref accum-array (1- end))))
             (pick (if (zerop max) (return-from select-random nil) (if fixed-pick fixed-pick (random max))))
             (found-total-index (array-search pick accum-array start end))
             (relative-index (if found-total-index (- found-total-index start) nil)))
        (print (list max pick))
        relative-index))))

;;; DEPRECATED, see +very-negative-prob-float+
;; (defmacro +small (arg1 arg2)
;;   `(handler-case (+ ,arg1 ,arg2)
;;      (floating-point-overflow () +most-negative-prob-float+)))
;; (defun +neg (arg1 arg2) ;sum of very small numbers to avoid underflows and signal errors that are time consuming
;;   (declare (optimize (speed 3) (safety 0)) (prob-float arg1 arg2))
;;   (if (and (< arg1 +very-negative-prob-float+) (< arg2 +very-negative-prob-float+))
;;       +most-negative-prob-float+
;;       (+ arg1 arg2)))


;;; functions to create the log versions of an array or a vector. Specialized for prob-float. Suitable to be a general utility
(defun log-vector (array)
  (declare ((prob-array (*)) array))
  (let ((out (make-array (array-dimensions array) :element-type 'prob-float)) (v +0-prob+))
    (declare ((prob-array (*)) out) (prob-float v))
    (dotimes (i (array-dimension array 0) out)
      (setf (aref out i) (progn (setq v (aref array i))
                                (if (zerop v) +very-negative-prob-float+ (the prob-float (log v))))))))
(defun log-array (array)
  (declare ((prob-array (* *)) array))
  (let ((out (make-array (array-dimensions array) :element-type 'prob-float)) (v +0-prob+))
    (declare ((prob-array (* *)) array) (prob-float v))
    (dotimes (i (array-dimension array 0) out)
      (dotimes (j (array-dimension array 1))
        (setf (aref out i j)
              (progn (setq v (aref array i j))
                     (if (zerop v) +very-negative-prob-float+ (the prob-float (log v)))))))))

;;; Copied from cl-utilities
(defun split-sequence (delimiter seq &key (count nil) (remove-empty-subseqs nil) (from-end nil) (start 0) (end nil) (test nil test-supplied) (test-not nil test-not-supplied) (key nil key-supplied))
  "Return a list of subsequences in seq delimited by delimiter.

If :remove-empty-subseqs is NIL, empty subsequences will be included
in the result; otherwise they will be discarded.  All other keywords
work analogously to those for CL:SUBSTITUTE.  In particular, the
behaviour of :from-end is possibly different from other versions of
this function; :from-end values of NIL and T are equivalent unless
:count is supplied. The second return value is an index suitable as an
argument to CL:SUBSEQ into the sequence indicating where processing
stopped."
  (let ((len (length seq))
        (other-keys (nconc (when test-supplied
                             (list :test test))
                           (when test-not-supplied
                             (list :test-not test-not))
                           (when key-supplied
                             (list :key key)))))
    (unless end (setq end len))
    (if from-end
        (loop for right = end then left
           for left = (max (or (apply #'position delimiter seq
                                      :end right
                                      :from-end t
                                      other-keys)
                               -1)
                           (1- start))
           unless (and (= right (1+ left))
                       remove-empty-subseqs) ; empty subseq we don't want
           if (and count (>= nr-elts count))
           ;; We can't take any more. Return now.
           return (values (nreverse subseqs) right)
           else
           collect (subseq seq (1+ left) right) into subseqs
           and sum 1 into nr-elts
           until (< left start)
           finally (return (values (nreverse subseqs) (1+ left))))
        (loop for left = start then (+ right 1)
           for right = (min (or (apply #'position delimiter seq
                                       :start left
                                       other-keys)
                                len)
                            end)
           unless (and (= right left)
                       remove-empty-subseqs) ; empty subseq we don't want
           if (and count (>= nr-elts count))
           ;; We can't take any more. Return now.
           return (values subseqs left)
           else
           collect (subseq seq left right) into subseqs
           and sum 1 into nr-elts
           until (>= right end)
           finally (return (values subseqs right))))))


(defun read-pair-observations-file (path)
  (with-open-file (stream path)
    (loop for line = (read-line stream nil)
       while line
       for (x_ y_ p) = (split-sequence #\Space line)
       for x = (mapcar #'parse-integer (split-sequence #\_ x_))
       for y = (mapcar #'parse-integer (split-sequence #\_ y_))
       with left-alphabet = (make-hash-table :test 'equalp)
       with right-alphabet = (make-hash-table :test 'equalp)
       do
         (loop for l in x do
              (setf (gethash l left-alphabet) nil))
         (loop for r in y do
              (setf (gethash r right-alphabet) nil))
       collecting (list (make-array (length x) :initial-contents x) (make-array (length y) :initial-contents y)) into observations
       finally (return (values observations (hash-table-size left-alphabet) (hash-table-size right-alphabet))))))
