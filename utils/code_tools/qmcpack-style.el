;;; qmpack-style.el -- defines a c-style for QMCPACK

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

;;  This attempts to give on the fly c-style that is as close as possible to the
;;  qmcpack coding standards and clang-format setup.
;;  You should still run clang-format on each file before pushing.

;;  I use this by
;;  (add-to-list 'load-path "/home/epd/qmcpack/utils/code_tools")
;;  (require 'qmcpack-style)
;;
;;  customize C Default Style if you want this style by default.
;;
;;  alternatively add a .dir-locals.el file to qmcpack/src
;;  by running this elisp and saving the .dir.locals.el buffer that appears
;;
;;  (let ((default-directory "~/qmcpack/src"))
;;     (add-dir-local-variable nil 'c-default-style "qmcpack"))
;;

;;; Code:
(defconst qmcpack-c-style
  '((c-basic-offset . 2)
    (c-offsets-alist
    (knr-argdecl-intro . +)
    (substatement-label . 2)
    (label . 2)
    (inline-open . 0)
    (func-decl-cont . +)
    (knr-argdecl . 0)
    (annotation-top-cont . 0)
    (annotation-var-cont . +)
    (member-init-cont . c-lineup-multi-inher)
    (inher-intro . +)
    (block-open . 0)
    (brace-list-open . 0)
    (brace-list-close . 0)
    (brace-list-intro . 0)
    (brace-list-entry . c-lineup-under-anchor)
    (brace-entry-open . 0)
    (statement-case-intro . +)
    (statement-case-open . 0)
    (case-label . 0)
    (do-while-closure . 0)
    (catch-clause . 0)
    (arglist-cont c-lineup-gcc-asm-reg 0)
    (cpp-macro-cont . +)
    (cpp-define-intro c-lineup-cpp-define +)
    (friend . 0)
    (objc-method-intro .
                     [0])
    (objc-method-args-cont . c-lineup-ObjC-method-args)
    (objc-method-call-cont c-lineup-ObjC-method-call-colons c-lineup-ObjC-method-call +)
    (extern-lang-open . 0)
    (module-open . 0)
    (composition-open . 0)
    (extern-lang-close . 0)
    (module-close . 0)
    (composition-close . 0)
    (inextern-lang . +)
    (inmodule . +)
    (incomposition . +)
    (template-args-cont c-lineup-template-args +)
    (inlambda . c-lineup-inexpr-block)
    (lambda-intro-cont . +)
    (inexpr-statement . +)
    (inexpr-class . +)
    (topmost-intro . 0)
    (namespace-open . 0)
    (innamespace . 0)
    (topmost-intro-cont . 0)
    (class-open . 0)
    (inclass . +)
    (defun-block-intro . +)
    (statement . 0)
    (inline-close . 0)
    (substatement . +)
    (substatement-open . 0)
    (statement-block-intro . +)
    (block-close . 0)
    (statement-cont . ++)
    (class-close . 0)
    (namespace-close . 0)
    (member-init-intro . ++)
    (defun-open . 0)
    (defun-close . 0)
    (stream-op . 3)
    (else-clause . 0)
    (arglist-intro . ++)
    (access-label . 0)
    (c . c-lineup-C-comments)
    (inher-cont . c-lineup-multi-inher)
    (string . -1000)
    (comment-intro . c-lineup-comment)
    (arglist-cont-nonempty . c-lineup-arglist)
    (arglist-close . c-lineup-close-paren)
    (cpp-macro . -1000)))
  "QMCPACK C/C++ programming Style.")
(c-add-style "qmcpack" qmcpack-c-style)

(provide 'qmcpack-style)
;;; qmcpack-style.el ends here
