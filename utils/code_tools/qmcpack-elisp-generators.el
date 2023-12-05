;;; qmpack-elisp-generators.el -- commands to generate C++ code

;;; License:
;; //////////////////////////////////////////////////////////////////////////////////////
;; // This file is distributed under the University of Illinois/NCSA Open Source License.
;; // See LICENSE file in top directory for details.
;; //
;; // Copyright (c) 2021 QMCPACK developers.
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

;;  To use you should load the functions
;;
;;  either put in your .emacs or
;;  run from *scratch* or by other means after emacs startup
;;
;;  I use this by
;;  (add-to-list 'load-path "/home/epd/qmcpack/utils/code_tools")
;;  (require 'qmcpack-elisp-generators)
;;
;;  then in your c++ mode buffer
;;  select the region with your variables
;;  then M-x qmcp-add-getters
;;  the generated code will appear after the end of the selected region in your buffer.
;;
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

(setq qmcp-variable-declaration-re "\\( *\\)\\([<>,A-Za-z:_&\\*]+ *[<>,A-Za-z:_&\\*]+\\)\\( +\\)\\([A-Za-z_0-9]+\\)_.*;")

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
              )
          (setq getters "")
          (while (re-search-forward qmcp-variable-declaration-re nil t)
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
          )
      (setq setters "")
      (while (re-search-forward qmcp-variable-declaration-re nil t)
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

(defun qmcp-convert-ptr-vector-ref-vector()
  "For std::vector<TYPE*> selected in REGION replace with RefVector<TYPE>."
  (interactive)
  (let ((leading-indent (qmcp-leading-indent)))
  (save-excursion
    (save-restriction
    (narrow-to-region (region-beginning) (region-end))
    (goto-char (point-min))
    (let ((my-type "")
          (post-fix "")
          (my-var "")
          (ref-vector-line "")
          (ptr-vector-declaration-re "std::vector<\\(.*?\\)\\**>\\([ &]*\\) *\\([A-Za-z0-9_]*\\)")
          )
      (setq these-vectors "")
      (while (re-search-forward ptr-vector-declaration-re nil t)
        (setq my-type (match-string 1))
        (setq post-fix (match-string 2))
        (setq my-var (match-string 3))
        (setq ref-vector-line (format "RefVector<%s>%s %s" my-type post-fix my-var))
        (replace-match ref-vector-line)))))))

(defun qmcp-revert-ptr-vector-ref-vector()
  "For std::vector<TYPE*> selected in REGION replace with RefVector<TYPE>."
  (interactive)
  (let ((leading-indent (qmcp-leading-indent)))
  (save-excursion
    (save-restriction
    (narrow-to-region (region-beginning) (region-end))
    (goto-char (point-min))
    (let ((my-type "")
          (post-fix "")
          (my-var "")
          (ref-vector-line "")
          (ptr-vector-declaration-re "RefVector<\\(.*?\\)>\\([ &]*\\) *\\([A-Za-z0-9_]*\\)")
          )
      (setq these-vectors "")
      (while (re-search-forward ptr-vector-declaration-re nil t)
        (setq my-type (match-string 1))
        (setq post-fix (match-string 2))
        (setq my-var (match-string 3))
        (setq ref-vector-line (format "std::vector<%s*>%s %s" my-type post-fix my-var))
        (replace-match ref-vector-line)))))))

(defun s-replace (old new s)
  "Replaces OLD with NEW in S."
  (declare (pure t) (side-effect-free t))
  (replace-regexp-in-string (regexp-quote old) new s t t))

(defun qmcp-add-setIfInInput()
  "For each C++ variable declaration in REGION write setIfInInput.
   The setIfInInput functions are written on starting on the line after the REGION and assume
   that the tag is the same as the var name less _."
  (interactive)
  (let ((leading-indent (qmcp-leading-indent)))
    (save-excursion
      (save-restriction
        (narrow-to-region (region-beginning) (region-end))
        (goto-char (point-min))
        (let ((setIIIs "")
              (indent "")
              (my-type "")
              (my-var "")
              (my-tag "")
              (setIII-line "")
              (first-line (line-beginning-position))
              )
          (while (re-search-forward qmcp-variable-declaration-re nil t)
            (message "match")
            (setq indent (match-string 1))
            (when (= first-line (line-beginning-position))
              (setq indent (concat indent leading-indent)))
            (setq my-type (match-string 2))
            (setq my-var (match-string 4))
            (setq my-tag (match-string 4))
            (s-replace "_$" "" my-tag)
            (setq setIII-line (format "setIfInInput(%s_, \"%s\");\n" my-var my-tag))
            (setq setIIIs (concat setIIIs setIII-line)))
          (goto-char (point-max))
          (insert "\n")
          (insert setIIIs))))))

(defun qmcp-add-enum-string-map()
  "Select all input enum classes into a single region to create string to value lookup map
   for use with assignAnyEnum. The  unordered_map is written on starting on the line
   after the REGION"
  (interactive)
  (let ((leading-indent (qmcp-leading-indent)))
    (save-excursion
      (save-restriction
        (narrow-to-region (region-beginning) (region-end))
        (goto-char (point-min))
        (let ((mapping-lines "")
              (elem-lines)
              (first-enum t)
              (line "")
              (indent "")
              (my-value "")
              (my-string "")
              (my-var "")
              (my-type "")
              (my-type-string "")
              (my-elems-begin (point-min))
              (end-of-enum (point-max))
              (mapping-line "")
              (first-line (line-beginning-position))
              (enum-class-declaration-re "\\( *\\)enum class \\([A-Z][A-Za-z1-9]+\\).*$")
              (enum-elem-declaration-re "\\( *\\)\\([A-Z][A-Za-z1-9_]+\\).*$"))
          (while (setq my-elems-begin (re-search-forward enum-class-declaration-re nil t))
            (setq indent (match-string 1))
            (setq my-type (match-string 2))
            (setq my-var (downcase (match-string 2)))
            (setq my-type-string (downcase my-type))
            (if first-enum
                (progn (setq mapping-line (format "%sinline static const std::unordered_map<std::string, std::any>\n%slookup_input_enum_value\n%s{\n" indent indent indent))
                 (setq mapping-lines (concat mapping-lines mapping-line))
                 (setq first-enum nil)) nil)
            (setq end-of-enum (re-search-forward "\\}; *$"))
            (goto-char my-elems-begin)
            (while (re-search-forward enum-elem-declaration-re end-of-enum t)
              (setq my-value (match-string 2))
              (setq my-string (downcase (match-string 2)))
              (setq line (format "%s  {\"%s-%s\", %s::%s}," indent my-type-string my-string my-type my-value))
              (push line elem-lines)))
          (setq line (pop elem-lines))
          (setq line (replace-regexp-in-string ",$" "" line nil nil))
          (push line elem-lines)
          (dolist (line (reverse elem-lines) mapping-lines)
            (setq mapping-lines (concat mapping-lines (format "%s\n" line))))
          (setq mapping-lines (concat mapping-lines (format "%s};\n" indent)))
          (goto-char (point-max))
          (insert "\n")
          (insert mapping-lines))))))

(provide 'qmcpack-elisp-generators)


;;; qmcpack-elisp-generators.el ends here
