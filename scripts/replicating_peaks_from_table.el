;; in-buffer creation of non-replicating peak tables

;; setup:
;; make variables for each table row
;; and a variable to dump table output into
(defun replicating-peak-table-generator ()
  "Generate replicating QTL peak tables in an org file."
  (interactive)
(let ((row-string)
      (table-list)
      (replicate-strings '("Replicate 1" "Replicate 2")))

  (dolist (r replicate-strings)
    (goto-char (point-min))
    (re-search-forward (concat r " Raw Table") nil t)
    (next-line)
    (while (org-at-table-p)
      (forward-char 1))
    ;; set up a heading to dump non-replicating peaks into
    (insert (concat "\n* " r " non-Replicating Peaks"))
    (goto-char (point-min))
    (re-search-forward (concat r " Final Table"))
    (replace-match (concat r " Replicating Peaks"))

  ;; go back to the table
    (goto-char (point-min))
    (re-search-forward (concat r " Raw Table") nil t)
    (while (not (org-at-table-p))
      (forward-char))
    (org-table-analyze)
    (org-table-goto-field "@2$2")

  ;; while in the table, check whether the peak replicated
  (while (org-at-table-p)
    (goto-char (line-beginning-position))
    (setq row-string (buffer-substring (point) (line-end-position)))
    (next-line 1)
    ;; if test
    ;; if the peak replicated, add the output to our table list
    ;; if not, do nothing
    (if (not
         (save-excursion
           (re-search-forward row-string (point-max) t)))
        (setq table-list (cons row-string table-list))
      (message "no"))
    )

  ;; reverse the list to get the peak order right
  (setq table-list (reverse table-list))

  ;;table creation
  ;; insert table header
  (re-search-forward
   (concat "\* " r " Replicating Peaks") (point-max))
  (insert "\n|----------+-----------+-----+--------+----------+------------+-----------+-------------|\n| reporter | replicate | chr |    LOD | delta_AF | left_Index | max_Index | right_Index |\n|----------+-----------+-----+--------+----------+------------+-----------+-------------|\n")
  (mapcar (lambda (arg) (insert (concat arg "\n"))) table-list)
  (setq table-list nil)
)
  )
)

;; export of peaks as .csv files
(defun all-peaks-export-from-org ()
  "Export replicating and non-replicating QTLs from an org file as a table."
  (interactive)
(let* ((reporter (substring (buffer-name (current-buffer)) 0 3))
       (replicate-strings '("Replicate 1" "Replicate 2"))
       (presence-strings '(" Replicating Peaks" " non-Replicating Peaks"))
       (full-strings))

  ;; set up string list
  (dolist (r replicate-strings)
    (dolist (p presence-strings)
      (setq full-strings (cons (concat r p) full-strings))))
  (setq full-strings (reverse full-strings))

  (dolist (f full-strings)
    (goto-char (point-min))
    (re-search-forward f nil)
    (while (not (org-at-table-p))
      (forward-char))
    (org-table-export
     (concat "~/data/illumina/2021.10.30_all_UPS_rdata/peaks/replicating_peak_tables/"
             reporter "_"
             (replace-regexp-in-string " " "_" f) ".csv"))
    ))
)
