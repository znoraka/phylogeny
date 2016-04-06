#lang racket

(struct distance-formula (global-func per-item-func) #:transparent)
(define d-f distance-formula)

(define (quantize n l)
  (define (q i)
    (floor (/ (+ i (floor (/ n 2))) n)))
  (define (iq i)
    (* n i))

  (map (λ (i)
         (add1 (iq (q i)))) l))

(define (random-values n)
  (let* ([l (build-list n (λ (i) (add1 (random 10000))))]
         [s (apply + l)])
    ;; (map (λ (i) (/ i s)) l)))
    l))

(define (distance dist l1 l2)
  (define (d a b)
    ((distance-formula-global-func dist) (apply + (map (λ (i j)
                                                         ((distance-formula-per-item-func dist) i j)) a b))))
  (if (= (length l1) (length l2))
      `(,(d l1 l2) ,(d l2 l1))
      '(+inf.0 +inf.0)))
  
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
  (let* ([max-value (apply max l)]
         [histo (make-vector (add1 max-value) 1)])
    (for-each (λ (i)
                (vector-set! histo i (add1 (vector-ref histo i)))) l)
    (vector->list histo)))

(define n 10000)

(define r1 (random-values n))
(define r2 (random-values 100))
(define r3 (random-values 100))

(define r4 `(65 ,(no-dup-q 65 (no-dup-q 85 r1))))
(define r5 `(65 ,(no-dup-q 65 r1)))
(define r6 `(35 ,(no-dup-q 35 r1)))
(define r7 `(85 ,(no-dup-q 85 r1)))
(define r8 `(65 ,(no-dup-q 65 (no-dup-q 75 (no-dup-q 85 r1)))))
(define r9 `(75 ,(no-dup-q 75 (no-dup-q 85 r1))))
(define r10 `(55 ,(no-dup-q 55 (no-dup-q 65 (no-dup-q 75 (no-dup-q 85 r1))))))

(define (can-be-child? child parent)
  (if (> (length (cdr child)) (length (cadr parent)))
      #f
      (let* ([q (min (car child) (car parent))]
             [parent (normalize (no-dup-q q (cadr parent)))]
             ;; [parent (normalize (cadr parent))]
             [child (normalize (cadr child))])
        (= 100
           (* 100.0
              (/ (count (λ (i)
                          (member i parent)) child)
                 (length child)))))))

(define (unit-test r1 r2 t1 t2)
  (define (from-string s)
    (eval (read (open-input-string s))))

  (define dists (distance kb
                          (cdadr r1)
                          (cdadr r2)))
  
  (define (f r1 r2 t d s)
    (let ([res (can-be-child? r1 r2)])
      (cond
        [(and res (not t)) (set! faux-positif (add1 faux-positif))]
        [(and (not res) t) (set! faux-negatif (add1 faux-negatif))]
        [(equal? res t) (set! correct (add1 correct))])
      ;; (displayln (~a s " => " (round d) " => " res "/" t " => " (if (equal? res t) #t "#f <-----")))))
      ))
  
  (f r1 r2 t1 (car dists) "r1 r2")
  (f r2 r1 t2 (cadr dists) "r2 r1"))

(define (make-test-data)
  (define r1 (random-values 10000))
  
  (define (make-compression-history n)
    (reverse (foldl (λ (i l)
                      (cons
                       (max 30 (- (car l) (* 5 (add1 (random 4)))))
                       l))
                    (list (- 100 (* 5 (random 4)))) (make-list n 0))))

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
    
    (unless (andmap (λ (i j)
                      (= i j)) (take h1 l) (take h2 l))
      (unit-test r1 r2 #f #f)))

  (if (= 0 (random 2))
      (test-false-data)
      (test-true-data)))

(define faux-positif 0)
(define faux-negatif 0)
(define correct 0)

(define (tests)
  (for* ([i (in-range 1000)])
    (make-test-data))

  (define total (+ faux-negatif faux-positif correct))
  
  (displayln (~a "faux positifs = " faux-positif " => " (* 100.0 (/ faux-positif total))))
  (displayln (~a "faux negatifs = " faux-negatif " => " (* 100.0 (/ faux-negatif total))))
  (displayln (~a "corrects = " correct " => " (* 100.0 (/ correct total)))))

  ;; (unit-test "r4" "r5" #f #f)
  ;; (unit-test "r4" "r7" #t #f)
  ;; (unit-test "r4" "r8" #f #f)
  ;; (unit-test "r8" "r9" #t #f)
  ;; (unit-test "r8" "r10" #f #t)
  ;; (unit-test "r10" "r9" #t #f))

