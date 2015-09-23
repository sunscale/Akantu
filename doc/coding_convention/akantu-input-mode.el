(defconst akantu-input-regex-section-opening "^[ \t]*\\([_A-Za-z]+\\)[ \t]+\\([_A-Za-z][_A-Za-z0-9]+\\)\\([ \t]+[_A-Za-z][_A-Za-z0-9]*\\)?[ \t]*\\[\\(#.*\\)?$")
(defconst akantu-input-regex-section-closing "^[ \t]*\\][ \t]*$")
(defconst akantu-input-regex-multiline "^.*\\\\[ \t]*\\(#.*\\)?$")
(defconst akantu-input-regex-variable-assign "^[ \t]*\\([A-Za-z][_A-Za-z0-9]*\\)[ \t]*=")

(defvar akantu-intput-mode-hook nil)
(defvar akantu-intput-mode-map
  (let ((map (make-keymap)))
    (define-key map "\C-j" 'newline-and-indent)
    map)
  "Keymap for Akantu input files major mode")

(defun akantu-input-indent-line ()
  "Indent current line as akantu input file."
  (interactive)
  (if (bobp)  ; Check for rule 1
      (indent-line-to 0)
    (let ((not-indented t) cur-indent)
      (setq cur-column nil)
      (setq cur-indent 0)
      (save-restriction
        (widen)
        (save-excursion
          (if (not (or (looking-back "^\\s-+") (bolp)))
              (progn
                (setq cur-column (+ (current-column) 1))
                (beginning-of-line)
                (search-forward-regexp "\\S-")
                (setq cur-column (- cur-column (current-column))))))
        (save-excursion
          (beginning-of-line)
          (forward-line -1)
          (if (looking-at akantu-input-regex-multiline)
              (progn
                (setq not_assignment t)
                (while not_assignment
                  (if (looking-at akantu-input-regex-variable-assign)
                      (setq not_assignment nil)
                    (progn
                      (forward-line -1)
                      (message "toto"))))
                (search-forward "=")
                (setq cur-indent (+ (+ cur-indent (current-column)) 1)))))
        (save-excursion
          (beginning-of-line)
          (if (looking-at akantu-input-regex-section-opening)
              (setq cur-indent (- cur-indent tab-width)))
          (while not-indented
            (if (looking-at akantu-input-regex-section-opening)
                (setq cur-indent (+ cur-indent tab-width))
              (if (looking-at akantu-input-regex-section-closing)
                  (setq cur-indent (- cur-indent tab-width))))
            (if (bobp)
                (setq not-indented nil))
            (forward-line -1))))
      (if (< cur-indent 0)
          (setq cur-indent 0))
      (indent-line-to cur-indent)
      (if cur-column
          (move-to-column (+ cur-column cur-indent))))))

(defconst akantu-input-font-lock-buffer
  (list
   '("\\<\\(global\\|material\\|model\\|mesh\\|heat\\|contact\\|friction\\|embedded_interface\\|rules\\|non_local\\|user\\)\\>"
     1 font-lock-keyword-face) ;; include definition
   '("^[ \t]*\\([A-Za-z][_A-Za-z0-9]*\\)[ \t]*=" 1 font-lock-variable-face)
   '("^[ \t]*\\([_A-Za-z]+\\)[ \t]+\\([_A-Za-z][_A-Za-z0-9]+\\)\\([ \t]+[_A-Za-z][_A-Za-z0-9]*\\)?[ \t]*\\[" 2 font-lock-reference-face)
   '("^[ \t]*\\([_A-Za-z]+\\)[ \t]+\\([_A-Za-z][_A-Za-z0-9]+\\)\\([ \t]+[_A-Za-z][_A-Za-z0-9]*\\)?[ \t]*\\[" 3 font-lock-type-face))
  "akantu-input mode syntax highlighting."
  )

(defvar akantu-input-mode-syntax-table nil "Syntax table for akantu-input-mode.")
(setq akantu-input-mode-syntax-table nil)

(defun akantu-input-mode ()
  "Major mode for editing Akantu input files"

  (interactive)
  (kill-all-local-variables)
  (setq major-mode 'akantu-input-mode)
  (setq mode-name "akantu-input")

  ;; Create the syntax table
  (setq akantu-input-mode-syntax-table (make-syntax-table))
  (set-syntax-table akantu-input-mode-syntax-table)
  (modify-syntax-entry ?_ "w" akantu-input-mode-syntax-table)
  (modify-syntax-entry ?# "<"  akantu-input-mode-syntax-table)
  (modify-syntax-entry ?\n ">" akantu-input-mode-syntax-table)

  ;; Setup font-lock mode, indentation and comment highlighting
  (set (make-local-variable 'indent-line-function) 'akantu-input-indent-line)

  (make-local-variable 'font-lock-defaults)
  (setq font-lock-defaults '(akantu-input-font-lock-buffer))

  (make-local-variable 'comment-start)
  (setq comment-start "#")
  (run-hooks 'akantu-input-mode-hook))

(provide 'akantu-input-mode)
