;;; qmpack-elisp-generators.el -- commands to generate C++ code

;;; License:
;; //////////////////////////////////////////////////////////////////////////////////////
;; // This file is distributed under the University of Illinois/NCSA Open Source License.
;; // See LICENSE file in top directory for details.
;; //
;; // Copyright (c) 2019 QMCPACK developers.
;; //
;; // File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
;; //
;; // File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
;; //////////////////////////////////////////////////////////////////////////////////////

;;; Commentary:

;;  This package provides a number of functions to generate C++ code from
;;  other C++ code selected in a region.  It assumes QMCPACK it is reading code
;;  following QMCPACK coding conventions and should produce code following
;;  the same conventions.

;;  There is more copy pasting than there should be feel, free to refactor

;;; Code:
(defun qmcp-leading-indent ()
  "Amount of first line outside of REGION.
when region begins without capturing the first lines indent grab it as
a correction string"
      (save-excursion
      (goto-char (region-beginning))
      (if (not(= (point) (line-beginning-position)))
          (make-string (- (point) (line-beginning-position)) ? )
        "")))

(defun qmcp-add-getters()
  "For each C++ variable declaration in REGION write getter.
The getter functions are written on starting on the line after the REGION."
  (interactive)
  (let ((leading-indent (qmcp-leading-indent)))
    (save-excursion
      (save-restriction
        (narrow-to-region (region-beginning) (region-end))
        (goto-char (point-min))
        (let ((getters "")
              (indent "")
              (my-type "")
              (my-var "")
              (getter-line "")
              (first-line (line-beginning-position))
              (variable-declaration-re "\\( *\\)\\([<>,A-Za-z:]+\\)\\( +\\)\\([A-Za-z_0-9]+\\)_.*;"))
          (setq getters "")
          (while (re-search-forward variable-declaration-re nil t)
            (message "match")
            (setq indent (match-string 1))
            (when (= first-line (line-beginning-position))
              (setq indent (concat indent leading-indent)))
            (setq my-type (match-string 2))
            (setq my-var (match-string 4))
            (setq getter-line (format "%s%s get_%s() const { return %s_; }\n" indent my-type my-var my-var))
            (setq getters (concat getters getter-line)))
          (goto-char (point-max))
          (insert "\n")
          (insert getters))))))

(defun qmcp-add-setters()
  "For each C++ variable declaration in REGION write setter.
The getter functions are written on starting on the line after the REGION."
  (interactive)
  (let ((leading-indent (qmcp-leading-indent)))
  (save-excursion
    (save-restriction
    (narrow-to-region (region-beginning) (region-end))
    (goto-char (point-min))
    (let ((setters "")
          (indent "")
          (my-type "")
          (my-var "")
          (setter-line "")
          (first-line (line-beginning-position))
          (variable-declaration-re "\\( *\\)\\([A-Za-z:]+\\)\\( +\\)\\([A-Za-z_0-9]+\\)_.*;"))
      (setq setters "")
      (while (re-search-forward variable-declaration-re nil t)
        (setq indent (match-string 1))
        (when (= first-line (line-beginning-position))
          (setq indent (concat indent leading-indent)))
        (setq my-type (match-string 2))
        (setq my-var (match-string 4))
        (setq setter-line (format "%svoid set_%s(%s %s) { %s_ = %s; }\n" indent my-var my-type my-var my-var my-var))
        (setq setters (concat setters setter-line)))
      (goto-char (point-max))
      (insert "\n")
      (insert setters))))))
(provide 'qmcpack-elisp-generators)


;;; qmcpack-elisp-generators.el ends here
