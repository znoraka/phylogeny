#lang racket

(require racket/file)
(require "scripts.rkt")

(struct distance-formula (global-func per-item-func) #:transparent)
(define d-f distance-formula)

(define (quantize n l)
  (define (q i)
    (floor (/ i n)))
    ;; (floor (/ (+ i (floor (/ n 2))) n)))
  (define (iq i)
    (* n i))

  (map (λ (i)
         ;; (values (q i))) l))
         (iq (q i))) l))

(define (random-values n)
  (let* ([l (build-list n (λ (i) (random 500)))]
         [s (apply + l)])
    ;; (map (λ (i) (/ i s)) l)))
    l))

(define (distance dist l1 l2)
  (define (d a b)
    ((distance-formula-global-func dist) (apply + (map (λ (i j)
                                                         ((distance-formula-per-item-func dist) i j)) a b))))

  (define l (min (length l1) (length l2)))
  ;; (if (= (length l1) (length l2))
  ;; `(,(d (take l1 l) (take l2 l)) ,(d (take l2 l) (take l1 l))))
  (d (take l1 l) (take l2 l)))
  ;;     '(+inf.0 +inf.0)))
  
(define euclid
  (d-f sqrt
       (λ (i j) (* (/ i j) (expt (- i j) 2)))))        

(define bhattacharyya
  (d-f log
       (λ (i j) (sqrt (* i j)))))

(define kb
  (d-f values
       (λ (i j) (* i (log (/ i j))))))

(define hellinger
  (d-f (λ (i) (sqrt (* 2 i)))
       (λ (i j) (expt (- (sqrt i) (sqrt j)) 2))))

(define hellinger-b
  (d-f (λ (i) (* 2 (sqrt (- 1 i))))
       (λ (i j) (sqrt (* i j)))))

(define s-hellinger
  (d-f (λ (i) (* 2 i))
       (λ (i j) (expt (- (sqrt i) (sqrt j)) 2))))

(define jeffreys
  (d-f values
       (λ (i j) (* (- i j) (log (/ i j))))))

(define kdiv
  (d-f values
       (λ (i j) (* i (log (/ (* 2 i) (+ i j)))))))
        
(define (no-dup-q q l)
  (filter-not zero? (remove-duplicates (quantize q l))))

(define (normalize l)
  (no-dup-q 100 l))

(define (make-histo l)
  (let* ([l (map abs l)]
         [max-value (apply max l)]
         [histo (make-vector (add1 max-value) 1)])
    (for-each (λ (i)
                (vector-set! histo i (add1 (vector-ref histo i)))) l)
    (vector->list histo)))

(define (make-distrib l)
  (map (λ (i)
         (/ i (length l))) l))

(define n 1000)
(define r1 (random-values n))

(define (can-be-child? child parent)
  (define (f c p q)
    (* 100.0
       (/ (count (λ (i)
                   (ormap (λ (pj)
                            (and (>= i (- pj (/ q 4)))
                                (<= i (+ pj (/ q 4))))) p))
                   ;; (member i p))
                   
                 c)
          (length c))))

  ;; (if (equal? child parent)
      ;; #f  
      ;; (if (> (length (cadr child)) (length (cadr parent)))
          ;; #f
          (let* ([q (min (car child) (car parent))]
                 ;; [parent (values (no-dup-q q (cadr parent)))]
                 ;; [parent1 (no-dup-q q (cadr parent))]
                 ;; [child1 (normalize (no-dup-q q (cadr child)))]
                 ;; [parent2 (normalize (no-dup-q q (cadr parent)))]
                 ;; [child2 (no-dup-q q (cadr child))]
                 [parent3 (cadr parent)]
                 [child3 (cadr child)]
                 ;; [d1 (distance kb child2 parent1)]
                 ;; [d2 (distance kb child2 parent2)]
                 [d3 (distance kb child3 parent3)]
                 ;; [v (f child1 parent1 q)]
                 ;; [v2 (f child2 parent2 q)]
                 ;; [v3 (f child3 parent3 q)]
                 ;; [v4 (f child2 parent1 q)]
                 )
            ;; (displayln (~a "v = " v3))
            ;; (displayln (~a v " " v2 " " v3 " " v4))
            ;; (displayln (~a v3 " " v2))
            ;; (displayln (~a  v2))
            ;; (displayln (~a d1 " " d2 " " d3))
            ;; (displayln d3)
            ;; (newline)
            ;; (>= v3 100)))
            ;; (>= v2 80)))
            ;; (and (< d1 0) (< d2 0))))
            (= d3 0)))
            ;; (> v3 98)))
            ;; (and (>= v 50) (>= v3 50) (>= v4 50))))
            ;; (and (>= v 50) (> v2 0))))))

(define (unit-test r1 r2 t1 t2)
  (define (from-string s)
    (eval (read (open-input-string s))))

  ;; (define dists (distance kb
  ;;                         (cdadr r1)
  ;;                         (cdadr r2)))

  ;; (define dists '(() ()))
  
  (define (f r1 r2 t d s)
    (let ([res (can-be-child? r1 r2)])
      (cond
        [(and res (not t)) (set! fpr (add1 fpr))]
        [(and (not res) t) (set! fnr (add1 fnr))]
        [(and (not res) (not t)) (set! tnr (add1 tnr))]
        [(and res t) (set! tpr (add1 tpr))])
      ))
  
  (f r1 r2 t1 0 "r1 r2")
  (f r2 r1 t2 0 "r2 r1"))

(define (make-test-data)
  ;; (define r1 (random-values 10000))
  
  (define (make-compression-history n)
    (foldl (λ (i l)
             (if (= (car l) 27)
                 l
                 (cons
                  ;;desendre d'au max 25 par q
                  ;; (max 30 (- (car l) (* 5 (add1 (random 4)))))
                  (min 27 (+ (car l) (add1 (random (- 27 (car l) )))))
                  l)))
           `(,(add1 (random 26)))
           (make-list n 0)))

  (define (make-data-list compression-history)
    (list (car (reverse compression-history))
          (foldl (λ (i l)
                   (no-dup-q i l)) r1 compression-history)))

  (define (test-true-data)
    ;;avoir entre 2 et 4 q
    (define h (make-compression-history (+ 1 (random 3))))
    (define r1 (make-data-list h))
    ;;prendre une partie de l'historique de compression
    (define r2 (make-data-list (take h (add1 (random (- (length h) 1))))))

    (unit-test r1 r2 #t #f))

  (define (test-false-data)
    (define h1 (make-compression-history (+ 1 (random 3))))
    (define h2 (make-compression-history (+ 1 (random 3))))
    (define l (min (length h1) (length h2)))
    (define r1 (make-data-list h1))
    (define r2 (make-data-list h2))

    ;;enlever les historiques identiques
    (unless (andmap (λ (i j)
                      (= i j)) (take h1 l) (take h2 l))
      (unit-test r1 r2 #f #f)))

  (if (= 0 (random 2))
      (test-false-data)
      (test-true-data)))

(define tpr 0)
(define tnr 0)
(define fnr 0)
(define fpr 0)

(define (tests)
  (for* ([i (in-range 1000)])
    (make-test-data))

  (define total (+ tpr tnr fnr fpr))
  
  (displayln (~a "faux positifs = " fpr " => " (* 100.0 (/ fpr total))))
  (displayln (~a "faux negatifs = " fnr " => " (* 100.0 (/ fnr total))))
  (displayln (~a "vrais positifs = " tpr " => " (* 100.0 (/ tpr total))))
  (displayln (~a "vrais negatifs = " tnr " => " (* 100.0 (/ tnr total)))))

(define (test-dcts-folder path)
  
  (define (data-from-image path)
    (define l (read (open-input-file (path->string path))))
    (cons (car l) `(,(make-distrib (make-histo (cadr l))))))
    ;; (cons (car l) `(,(map abs (filter-not zero? (remove-duplicates (cadr l)))))))
    ;; (cons (car l) `(,(map (λ (i) (if (zero? i) 1 (abs i))) (cadr l)))))
    
  (define (file-name path)
    (car (string-split (car (reverse (string-split (path->string path) "/"))) ".")))
  
  (define files (cdr (find-files (λ (i) (= 1 1)) path)))

  (for-each (λ (i)
              (displayln (~a "parents of " (file-name i) " "
                             (foldl (λ (j l)
                                      (if (equal? i j)
                                          l
                                          (let ([image-i (data-from-image i)]
                                                [image-j (data-from-image j)])
                                            (display (~a `(,(file-name i) " - " ,(file-name j)) " "))
                                            (if (can-be-child? image-i image-j)
                                                (cons (file-name j) l) 
                                                l))))
                                    ' () files)))) files))

(define test-folder "/home/noe/Documents/dev/phylogeny/tree-generator/data/sameQ/dcts")

(define (test-dct-quantizations quantization-steps path)
  (define (make-data-list compression-history dct-coeffs)
    (list (car (reverse compression-history))
          (foldl (λ (i l)
                   (quantize i l)) dct-coeffs quantization-steps)))

  (define (data-from-image path)
    (cadr (read (open-input-file path))))

  (make-data-list quantization-steps (data-from-image path)))

(define (diff-between-dct-coeffs path1 path2)
  (define (data-from-image path)
    (cadr (read (open-input-file path))))

  (map (λ (i j)
         (- i j)) (data-from-image path1) (data-from-image path2)))


;; estimates the parents of the jpeg images set in path
(define (estimate-parents path)
  (define files (map path->string (find-files (λ (i)
                                                (regexp-match "^.*\\.jpg$" (path->string i))) path)))

  (define (file-name path)
    (car (string-split (car (reverse (string-split path "/"))) ".")))


  (for-each (λ (i)
              (displayln (~a "parents of " (file-name i) " "
                             (foldl (λ (j l)
                                      (if (equal? i j)
                                          l
                                          (let* ([quality-i (get-image-quality i)]
                                                 [quality-j (get-image-quality j)])
                                            (if (< quality-i quality-j)
                                                (let* ([image-i `(0 ,(make-distrib (make-histo (cadr (get-dct-coefficients i)))))]
                                                       [a (compress j "out.jpg" quality-i)]
                                                       [image-j `(0 ,(make-distrib (make-histo (cadr (get-dct-coefficients "out.jpg")))))])
                                                  ;; (display (~a `(,(file-name i) " - " ,(file-name j)) " "))
                                                  (if (can-be-child? image-i image-j)
                                                      (cons (file-name j) l) 
                                                      l))
                                                l)))) '() files)))) files))
