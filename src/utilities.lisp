;; Author: Juan Miguel Cejuela
;; Created: Wed Jul  9 19:13:54 2008 (CEST)
;; Last-Updated: 2011-08-11
;;           By: Juan Miguel Cejuela
;;     Update #: 15

(in-package :cl-hmm)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;fastest search in a prob-float array, (see jmc.cl.utils)
(def-lin-search array-search prob-float (simple-array prob-float) >=)

;;; DEPRECATED, see +very-negative-prob-float+
;; (defmacro +small (arg1 arg2)
;;   `(handler-case (+ ,arg1 ,arg2)
;;      (floating-point-overflow () +most-negative-prob-float+)))
;; (defun +neg (arg1 arg2) ;sum of very small numbers to avoid underflows and signal errors that are time consuming
;;   (declare (optimize (speed 3) (safety 0)) (prob-float arg1 arg2))
;;   (if (and (< arg1 +very-negative-prob-float+) (< arg2 +very-negative-prob-float+))
;;       +most-negative-prob-float+
;;       (+ arg1 arg2)))


;;; functions to create the log versions of an array or a vector. Specialized for prob-float. Suitable to be an utility
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
