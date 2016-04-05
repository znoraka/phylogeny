#lang racket

(define (quantize n l)
  (define (q i)
    (floor (/ (+ i (floor (/ n 2))) n)))
  (define (iq i)
    (* n i))

  (map (λ (i)
         (iq (q i))) l))

(define (test n1 n2 l)
  (let ([q (quantize n2 (quantize n1 l))])
    (displayln q)
    (map (λ (i)
           (modulo i n1)) q)))

(define (random-values n)
  (let* ([l (build-list n (λ (i) (add1 (random 100))))]
         [s (apply + l)])
    (map (λ (i) (/ i s)) l)))
;; l))

(struct distance-formula (global-func per-item-func) #:transparent)
(define d-f distance-formula)

(define (distance dist l1 l2)
  (define (d a b)
    ((distance-formula-global-func dist) (apply + (map (λ (i j)
                                                         ((distance-formula-per-item-func dist) i j)) a b))))

  (if (= (length l1) (length l2))
      `(,(d l1 l2) ,(d l2 l1))
      `(,(/ 0.0 0.0) ,(/ 0.0 0.0))))
  

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
        
(define r1 (random-values 100))
(define r2 (random-values 100))
(define r3 (random-values 100))

(define (no-dup-q q l)
  (remove-duplicates (quantize q l)))

(define (normalize l)
  (no-dup-q 100 l))

(define r4 `(65 ,(no-dup-q 65 (no-dup-q 85 (range 10000)))))
(define r5 `(65 ,(no-dup-q 65 (range 10000))))
(define r6 `(35 ,(no-dup-q 35 (range 10000))))
(define r7 `(85 ,(no-dup-q 85 (range 10000))))
(define r8 `(65 ,(no-dup-q 65 (no-dup-q 75 (no-dup-q 85 (range 10000))))))
(define r9 `(75 ,(no-dup-q 75 (no-dup-q 85 (range 10000)))))
(define r10 `(55 ,(no-dup-q 55 (no-dup-q 65 (no-dup-q 75 (no-dup-q 85 (range 10000)))))))

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
                          (cdadr (from-string r1))
                          (cdadr (from-string r2))))
  
  (define (f r1 r2 t d)
    (let ([res (from-string (~a `(can-be-child? ,r1 ,r2)))])
      (displayln (~a r1 " " r2 " => " (round d) " => " res "/" t " => " (if (equal? res t) #t "#f <-----")))))
  (f r1 r2 t1 (car dists))
  (f r2 r1 t2 (cadr dists)))

(define (tests)
  (unit-test "r4" "r5" #f #f)
  (unit-test "r4" "r7" #t #f)
  (unit-test "r4" "r8" #f #f)
  (unit-test "r8" "r9" #t #f)
  (unit-test "r8" "r10" #f #t)
  (unit-test "r10" "r9" #t #f))

