#lang racket

(provide get-image-quality)
(provide compress)
(provide get-dct-coefficients)

(define (exec-script script)
  (define tmp-port (current-output-port))
  (define op1 (open-output-string))
  (current-output-port op1)
  (system script)
  (define out (string-replace (get-output-string op1) "\n" ""))
  (current-output-port tmp-port)
  out)

(define (get-image-quality path)
  (string->number (exec-script (~a "identify -verbose " path " | grep \"Quality\" | cut -d \" \" -f4"))))

(define (compress input-path output-path quality)
  (exec-script (~a "convert -quality " quality " " input-path " " output-path)))

(define (get-dct-coefficients path)
  (read (open-input-string (exec-script (~a "../c++/bin/dctextractor " path)))))
