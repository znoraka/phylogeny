#lang racket

(define (parse-file path)
  (define (extract-number s)
    (last (string-split s " ")))
  
  (define file (open-input-file path))

  (map string->number (build-list 7 (λ (n) (extract-number (read-line file))))))

(define (compute-results path)
  (define mean-error '())
  (define overestimate '())
  (define underestimate '())
  (define roots '())
  (define edges '())
  (define leaves '())
  (define ancestry '())

  
  )