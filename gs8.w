\documentclass[fleqn]{article}

\topmargin=-0.3in
\textheight=8in
\oddsidemargin=0in
\textwidth=6.5in

\usepackage{indentfirst}
\usepackage{framed}
\usepackage[usenames,dvipsnames]{color}
\usepackage{makeidx}
\usepackage{asymptote}

\definecolor{darkgreen}{rgb}{0.0,0.5,0.0}
\definecolor{darkorange}{rgb}{0.5,0.5,0.0}
\definecolor{magenta}{rgb}{0.5,0.0,0.5}

\makeatletter
\newcommand{\verbatimfont}[1]{\def\verbatim@font{#1}}%
\makeatother

\setlength{\parindent}{0.6cm}

\newcommand{\refn}[1]{\raisebox{1ex}{\footnotesize{#1}}}
\makeindex

\begin{document}

\newpage

\setcounter{page}{1}
\renewcommand{\baselinestretch}{1.4142} \small\normalsize
\section*{Simulation generator}

\section*{Implementation}
The language used for this version of the program is Common LISP.  

\section*{Utility functions and random number generation}

Initialization:
<<>>=
(setf *read-default-float-format* 'double-float)
@

Debug function:
<<>>=
(defun shove (m x) (format t "~A: ~A~%" m x) (finish-output NIL))
(defun sho (x) (format t "arg=~A~%" x) (finish-output NIL) x)
@

\subsection*{Normalization}
The function {\tt nml} accepts a list of numbers, and returns a list of the numbers normalized by
dividing by their sum.\index{nml}
<<>>=
(defun nml (aa) 
  (let ((asum (reduce '+ aa)))
    (if (equal asum 0) (error "zerodivide in nml"))
    (mapcar #'(lambda (u) (/ u asum)) aa)))
@

\subsection*{Lists}
Look for duplicates, based on {\tt eq}:
<<>>=
(defun any-dup-p (alist)
  (let ((ht (make-hash-table :test #'eq))
        (ans NIL))
    (tagbody
      (dolist (ii alist)
        (if (null (gethash ii ht))
            (setf (gethash ii ht) 1)
            (progn (setf ans T) (go quit-looking))))
      quit-looking
    )
    ans)) 
@

<<>>=
(defun vector-to-list (v)  ; refactor to eliminate this
  (let ((alist ()))
    (dotimes (ii (length v))
      (setf alist (cons (aref v ii) alist)))
    alist))
@

\subsection*{Park-Miller RNG}
A simple Park-Miller linear congruential generator is used for the time being.

The generator is initialized by {\tt setseed}.  The main function is {\tt nextrand}, which advances the generator.  Below,
an arbitrary value for the seed is used to initialize it.\index{setseed}\index{nextrand}
<<>>=
(let ((curseed 0))
  (defun setseed (newseed) 
    (if (<= newseed 0) 
        (error "seed must be positive"))
    (setf curseed newseed)
    T)
  (defun reportseed () curseed)
  (defun nextrand ()
    (setf curseed (mod (* 16807 curseed) 2147483647))
    curseed))
(setseed 6879098)
@

The function {\tt runif} called with no arguments yields a uniform variate between 0 and 1.\index{runif}
<<>>=
(defun runif () 
  (let ((ans (/ (nextrand) 2147483647.0)))
    ans))
(defun runif-range (low upp)
  (+ low (* (runif) (- upp low))))
@

The function {\tt rexp} with a single argument for the rate yields a single exponential
random variate with the given rate.\index{rexp}
<<>>=
(defun rexp (the-rate)
  (let ((uu (runif))
        (ans NIL))
    (cond ((<= uu 0)
           (error "zero from (runif)"))
          ((< (abs the-rate) 1e-20)
           (setf ans 1e20))  ; arbitrary big number
          (t 
           (setf ans (/ (- (log uu)) the-rate))))
    ans))
@

The function {\tt rbern-p} takes a probability and returns a Boolean which is {\tt T} with
the given probability and {\tt NIL} otherwise; no range checks are done.\index{rbern-p}
<<>>=
(defun rbern-p (prob) 
  (let ((ans (< (runif) prob)))
    ans))
@

The function {\tt rint} generates a discrete uniform
random integer between 0 and {\tt nn}-1.  The upper
endpoint is NOT a possible value.
<<>>=
(defun rint (nn)
  (let* ((tmp (reportseed))
         (urand (runif))
         (ans (floor (* nn urand))))
    (if (>= ans nn) 
        (rint nn)    
        ans)))
@

The function {\tt rint-range} generates a discrete uniform
from {\tt :lower} to {\tt :upper} {\it inclusive}.
<<>>=
(defun rint-range (&key lower upper)
  (if (< upper lower) (error "bad call to rint-range"))
  (let ((ns (+ (- upper lower) 1)))
    (+ lower (rint ns))))
@

The function {\tt rmultinom} takes a list of probabilities, and returns a single multinomial draw from
them: the probability of returning a {\tt 0} is given by the value of the first element, etc.  The list
is normalized by the function, input values out of range signal error. \index{probcheck}\index{rmultinom} 
<<>>=
(defun probcheck (plist)
  (cond ((> (apply #'max plist) 1) (error "mishap: prob > 1"))
        ((< (apply #'min plist) 0) (error "mishap: prob < 0"))
        (t plist)))
(defun rmultinom-aux (uu nn cumul alist)
  (cond ((null alist) nn)
        (t (let ((new (+ (car alist) cumul)))
             (cond ((< uu new) nn)
                   (t (rmultinom-aux uu (+ nn 1) new (cdr alist))))))))
(defun rmultinom (problist0)
  (let ((problist (nml (probcheck problist0)))
        (uu (runif))
        (ans NIL))
    (cond ((< uu (car problist)) 
           (setf ans 0))
          (t 
           (setf ans (rmultinom-aux uu 1 (car problist) (cdr problist)))))
    ans))
@

This selects a random member of a list.
<<>>=
(defun rselect (alist)
  (if (null alist) (error "null list in alist"))
  (let ((ans (nth (rint (length alist)) alist)))      ; %%%^^^
    ans))
(defun rselect-l (alist len)
  (nth (rint len) alist))
@

This function omits an element from a list.\index{omit}
<<>>=
(defun omit (n alist)
  (cond ((< n 0) (error "negative index"))
        ((equal n 0) (cdr alist))
        ((>= n (length alist))                        ; %%%^^^
         (error "out of range in omit"))
        (t (append (subseq alist 0 n ) (nthcdr (+ 1 n) alist)))))  ; from stackoverflow
@

This takes a sample without replacement of a list.\index{rsample} %%%!!! untested!
<<>>=
(defun rsample-aux (alist nsamp curlen cumul)
  (cond ((<= nsamp 0) cumul)
        (t (let* ((random-integer (rint curlen))
                  (picked (nth random-integer alist))
                  (newlist (omit random-integer alist)))
             (rsample-aux newlist (- nsamp 1) (- curlen 1) (cons picked cumul))))))
(defun rsample (alist nsamp)
  (let ((curlen (length alist)))                      ; %%%^^^
    (cond ((> nsamp curlen) (error "list too short"))
          (t (rsample-aux alist nsamp curlen ())))))
@

This takes a sample without replacement, based on weights.  It may be the
Wallenius noncentral hypergeometric.  Light objects weight 1; heavy objects weigh {\tt wt}.
<<>>=
(defun rsample-weighted (nsamp heavy light wt)
  (rsample-weighted-aux nsamp (shuffle heavy) (shuffle light) (length heavy) (length light) wt ()))
(defun rsample-weighted-aux (nsamp heavy light nh nl wt cumul)
  (cond ((and (> nsamp 0) (equal nh 0) (equal nl 0))
         (error "ran out of things to sample in rsample-weighted-aux"))
        ((<= nsamp 0)
         cumul)
        (t (let* ((heavy-list-weight (* wt nh))
                  (ph (/ heavy-list-weight (+ heavy-list-weight nl))))
             (if (or (< ph 0) (> ph 1))
                 (error "bad probability of heavy list"))
             (if (rbern-p ph)
                 (rsample-weighted-aux 
                  (- nsamp 1)
                  (cdr heavy)
                  light
                  (- nh 1)
                  nl 
                  wt 
                  (cons (car heavy) cumul)) 
                 (rsample-weighted-aux
                  (- nsamp 1)
                  heavy
                  (cdr light)
                  nh
                  (- nl 1)
                  wt
                  (cons (car light) cumul)))))))
@

This next function generates a binomial sample of the input list, given
the probability {\tt prob}.  The order is not preserved.\index{binomsamp}
<<>>=
(defun binomsamp (alist prob)
  (let ((cumul ()))
    (dolist (cur alist)
      (if (rbern-p prob) 
          (push cur cumul)))
    cumul))
@

This function shuffles a list.\index{shuffle}
It is redundant with rsample using the whole length of the list.
<<>>=
(defun shuffle-aux (ans len alist &optional (dbg NIL))
  (if (null alist) 
      ans
    (let* ((pos (rint len)))
      (if (>= pos len)
          (progn
           (error "shuffle-aux")))
      (if (and dbg (not (equal len (length alist))))
          (progn
           (error "shuffle-aux")))
      (shuffle-aux (cons (nth pos alist) ans) (- len 1) (omit pos alist)))))
(defun shuffle (xlist)
  (let ((ans (shuffle-aux () (length xlist) xlist))) 
    ans))
@

This shuffles a vector destructively:
<<>>=
(defun destructive-vector-shuffle (vec len) 
  (let ((rp NIL)
        (tmp NIL))
    (do ((cur (- len 1) (- cur 1)))
        ((<= cur 0) t)
      (setf rp (rint (+ 1 cur)))
      (setf tmp (aref vec cur))
      (setf (aref vec cur) (aref vec rp))
      (setf (aref vec rp) tmp))
    vec)) 
@

This just creates a sequence from 0, of the given length.
<<>>=
(defun vseq (nlen) 
  (let ((newa (make-array (list nlen))))
    (dotimes (ii nlen)
      (setf (aref newa ii) ii))
    newa))
@

\subsection*{Unique number generators}
It is helpful to have a unique number associated with every IDU agent, needle, and location (separately).
The functions {\tt new-idu-id}, {\tt new-needle-id}, and {\tt new-location-id} called with no arguments
yield such numbers.\index{id-maker}
<<>>=
(defun id-maker ()
  (let ((cur 0))
    (lambda ()
      (incf cur))))
; in poplog clisp the function print forms #<FUNCTION ...> aren't right but 
;   seem to work ok
(setf (symbol-function 'new-idu-id) (id-maker))
(setf (symbol-function 'new-needle-id) (id-maker))
(setf (symbol-function 'new-location-id) (id-maker))
@

\subsection*{Other useful functions}
We do not need both these:
<<>>=
(defun rep (obj ntimes)
  (let ((ans ()))
    (dotimes (ii ntimes)
      (setf ans (cons obj ans)))
    ans))
(defun dup (a x)
  ; duplicate a x times
  (make-sequence 'list x :initial-element a))
@

The function {\tt seq1} generates a sequence from {\tt start} to {\tt end} INCLUSIVE of the endpoints, in steps
of 1.\index{iseq}
<<>>=
(defun seq1 (start end)
  (cond ((>= start end) 
         (error "wrong way in seq"))
        (t (let ((ans ()))
             (dotimes (ii (+ (- end start) 1))
               (setf ans (cons (+ ii start) ans)))
             (reverse ans)))))
@

The following function generates sequences of integers; \verb:a:
and \verb:b: are the beginning and end of the sequence; the answer
is a list.  This should be the same as {\tt iseq}.\index{seq}
<<>>=
(defun seq (a b)
  (cond ((< a 0) (error "Negative first argument in seq."))
        ((< b 0) (error "Negative second argument in seq."))
        ((> a b) (error "First argument greater than second in seq."))
        (t (do* ((cur b (1- cur))
                 (ans (list b) (cons cur ans)))
             ((<= cur a) ans)))))
@

<<>>=
(defun fseq (start by len)
  (if (< (abs by) 1e-8)
      (error "bad by in fseq"))
  (loop for ii in (seq 0 (- len 1))
    collect (+ start (* ii by))))
@

Concatenate list (replace by builtin):\index{concatlist}
<<>>=
(defun concatlist-aux (inlist accumans sep)
  (cond ((null inlist) accumans)
        ((null (cdr inlist)) (concatenate 'string accumans sep (car inlist)))
        (t 
         (concatlist-aux (cdr inlist) 
                         (concatenate 'string accumans sep (car inlist)) sep))))
(defun concatlist (inlist sep)
  (concatlist-aux (mapcar #'(lambda (x) (format nil "~A" x)) inlist) "" sep))
@

Pull the first {\tt nn} objects from {\tt alist} satisfying {\tt pred}; return the objects 
and the remaining part of the list.\index{retrieve-group}
<<>>=
(defun retrieve-group-aux (ok nogood alist nn pred)
  (cond ((and (null alist) (> nn 0))
         (values ok (append (reverse nogood) alist) NIL))
        ((equal nn 0)
         (values ok (append (reverse nogood) alist) T))
        ((funcall pred (car alist))
         (retrieve-group-aux (cons (car alist) ok) nogood (cdr alist) (- nn 1) pred))
        (t
         (retrieve-group-aux ok (cons (car alist) nogood) (cdr alist) nn pred))))
(defun retrieve-group (alist nn pred)
  (cond ((< nn 0) (error "negative nn"))
        (t (retrieve-group-aux () () alist nn pred))))
@

\newpage
\section*{Code generation}
If I were to rewrite this, the code would be simpler and would depend on 
database methods.

\subsection{More utilities}
Both \verb+keys+ and \verb+target+ must be sorted.  Function returns true
if all elements of \verb+keys+ are in \verb+target+.
<<>>=
(defun list-includes-p (keys targets)
  (do ((curkeys keys (if (equalp (car curkeys) (car curtargets))
                         (cdr curkeys)
                       curkeys))    
       (curtargets targets (cdr curtargets)))
    ((or (null curkeys) (null curtargets)) (null curkeys))))
@

Determines whether two lists are equal, using \verb+equal+:
<<>>=
(defun equal-lists? (alist1 alist2)
  (do ((cur1 alist1 (cdr cur1))
       (cur2 alist2 (cdr cur2))
       (ans T (and ans (equal (car cur1) (car cur2)))))
    ((or (null cur1) (null ans)) ans)))
@

<<>>=
(defun list-of-symbols-p (xx)
  (if (consp xx)
      (do ((cur xx (cdr cur))
           (ans (if (null xx) NIL T) (and ans (symbolp (car cur)))))
        ((or (null cur) (null ans)) ans))
    NIL))
@

This produces a pair of lists from {\tt alist}, based on the predicate function {\tt pred}.
<<>>=
(defun split-list (pred alist)
  (cond ((null alist) (list () ()))
        (t (do* ((top      (first alist)      
                           (first cur))
                 (predcur? (funcall pred top) 
                           (funcall pred top))
                 (yesses   (cond (predcur? (list top))
                                 (t ()))
                           (cond (predcur? (cons top yesses))
                                 (t yesses)))
                 (nos      (cond (predcur? ())
                                 (t (list top)))
                           (cond (predcur? nos)
                                 (t (cons top nos))))
                 (cur      (rest alist) 
                           (rest cur)))
               ((null cur) (list (reverse yesses)
                                 (reverse nos)))))))
@

Functional ``and'':
<<>>=
(defun f/and (blist)
  (do ((cur blist (cdr cur))
       (ans (if (null blist) NIL T) (and ans (car cur))))
    ((or (null cur) (null ans)) ans)))
(defun f/or (blist)
  (cond ((null blist) NIL)
        ((car blist) T)
        (t (f/or (cdr blist)))))
@

<<>>=
(defun any (afn alist)
  (find-if #'(lambda (x) (funcall afn x)) alist))
(defun all (afn alist)
  (not (find-if-not afn alist)))
@

\subsection*{Names}
In this section, primitives for manipulating object names are developed.

For this version of the simulation generator, object names involve 
symbols.  Each symbol may be described as follows:
\begin{verbatim}
id.string ::= [a-z][a-z0-9]*
id ::= id.string | id.string\[[0-9][0-9]*\]
gname ::= id | id-gname
\end{verbatim}
For example, \verb+aa+, \verb+a4+, \verb+bb[5]+, \verb+aa-zz[4]+, 
\verb+a9-vv[6]-nn-z7+ are all allowable names.  Notionally, these
represent sets; \verb+aa-bb+ and \verb+bb-aa+ represent the
same thing. In fact \verb+bb-aa+ would not occur; internally these
names are always sorted.

The function \verb+cat+ is needed for printing; \verb+separator+ is a 
string, usually a dash or underscore; \verb+args+ is either a symbol,
a string, or a list of symbols and strings.  The output is a string.  
<<>>=
(defun cat (args)
  (do ((ans "" (concatenate 'string ans (symbol-name (car cur))))
       (cur (if (consp args) args (list args)) (cdr cur)))
    ((null cur) ans)))
@

The function \verb+show+ takes either a \verb+gname+, a \verb+transition+,
a \verb+definition+, a string or symbol, or a list of these, and returns
a string representing the input.
<<>>=
(defun show (x)
  (cond ((consp x) (concatenate 'string
                                  (cat "" x)
                                  " "))
        ((definition-p x) (show (definition-name x)))
        ((transition-p x) (show (transition-name x)))
        ((atom x) (concatenate 'string (string x) " "))
        (t (error "can't show"))))
@

First, the constructor \verb+new-gname+ accepts a list of symbols:
<<>>=
(defun new-gname (xx)
  (let* ((thelist (cond ((symbolp xx) (list xx))
                        ((stringp xx) (list (intern xx)))
                        ((consp xx)
                         (sort (copy-list xx) #'string<))
                        (t 
                         (prin1 xx) 
                         (error "new-gname: argument not list of symbols"))))
         (newsym (intern (cat thelist))))
    (setf (get newsym 'namelist) thelist)
    newsym))
@

\subsection*{Special syntax}

The next functions are modified from Slade.  They allow us to use
the syntax \verb:$xx: to refer to the state \verb+xx+.
<<>>=
(defvar gs-readtable
  (copy-readtable *readtable*))
(defun new-readtable (new)
  (setf *old-readtable* *readtable*)
  (setf *readtable* new))
(defun old-readtable ()
  (setf *readtable* *old-readtable*))
(set-macro-character #\$
  #'(lambda (stream char)
      (declare (ignore char))
      (list 'quote (new-gname (read stream t nil t))))
  t gs-readtable)
(new-readtable gs-readtable)
@

The next merge function does not change either of the original names,
but creates a new name.  The function \verb+merge-names-fn+ includes
a second return value, which is \verb+T+; this function is used in
the cross-stratification functions.
<<>>=
(defun merge-names (objname1 objname2)
  (cond ((or (eq objname1 '_SINK)
             (eq objname1 '_SOURCE)) objname1)
        ((or (eq objname2 '_SINK)
             (eq objname2 '_SOURCE)) objname2)
        (t (new-gname (merge 'list (copy-list (get objname1 'namelist))
                                   (copy-list (get objname2 'namelist))
                                   #'string<)))))
(defun merge-names-fn (objname1 objname2)
  (values (merge-names objname1 objname2) T))
(defun merge-names-excluding (objname1 objname2 excludefn)
  (cond ((or (eq objname1 '_SINK)
             (eq objname1 '_SOURCE)) objname1)
        ((or (eq objname2 '_SINK)
             (eq objname2 '_SOURCE)) objname2)
        (t (let* ((merged (merge 'list (copy-list (get objname1 'namelist))
                                       (copy-list (get objname2 'namelist))
                                       #'string<))
                  (newlist (remove-if excludefn merged)))
             (new-gname newlist)))))
@

<<>>=             
(defun j~ (objname1 alist)
  (if (null alist) 
      objname1
    (merge-names objname1 (j~ (car alist) (cdr alist)))))
(defun ~ (objname1 &rest args)
  (j~ objname1 args))
(defun source-p (oname)
  (find "_SOURCE" (get oname 'namelist) :test #'string=))
(defun sink-p (oname)
  (find "_SINK" (get oname 'namelist) :test #'string=))
(defun special-p (oname)
  (or (find "_SOURCE" (get oname 'namelist) :test #'string=)
      (find "_SINK" (get oname 'namelist) :test #'string=)))
@

Determine whether all the strings in \verb+oname1+ are also in \verb+oname2+:
<<>>=
(defun name-included-p (oname1 oname2)
  (list-includes-p (get oname1 'namelist) (get oname2 'namelist)))
@

Determine whether any of the names in \verb+oname1+ are included in
\verb+oname2+:
<<>>=
(defun any-included-p (onames1 oname2)
  (find-if #'(lambda (x) (name-included-p x oname2)) onames1))
@

The function \verb+same-name-p+ returns true if the arguments are
the same name, and raises an error if the arguments are not names.
<<>>=
(defun same-name-p (oname1 oname2)
  (eq oname1 oname2))
@

We could make the separator something other than a dash, but we don't:
<<>>=
(defun to-string (aname)
  (cat aname))
@

A vector quantity may be referred to as \verb+vv[a..b]+; such a
definition should always be included before any use of the vector or
any elements of it.  Once defined, a vector may be referred to 
by \verb+vv[]+ for example.  
To refer to a single vector element, you may write for instance
\verb+aa[3]+ to refer to element 3 of the vector \verb+aa+.  When the
index comes from evaluating an expression, the function \verb+@+ is
useful.  The input to \verb+@+ is a string and a number; the output is
a \verb+gname+:
<<>>=
(defun @ (anamestring arg)
  (new-gname
   (list 
    (intern 
     (map 'string #'char-upcase
      (remove (coerce " " 'character)
              (cond ((consp arg) 
                     (@v (format nil "~A[~A" anamestring (car arg)) (cdr arg)))
                    (t (format nil "~A[~A]" anamestring arg)))))))))
@

The function \verb+@v+ is a utility function used by \verb+@+.  It takes
a string \verb+cname+ and a list \verb+arg+; the result is a string which
consists of the \verb+cname+, and then \verb+,x+ repeated for each \verb+x+
in the \verb+arg+ list.  For example:
\begin{verbatim}
> (@v 'aa '(4 3))
"AA,4,3]"
> (@v "aa[1" '(4 3))
"aa[1,4,3]"
\end{verbatim}
Here is the function itself:
<<>>=
(defun @v (cname arg)
  (cond ((null arg) (format nil "~A]" cname))
        (t (@v (format nil "~A,~A" cname (car arg)) (cdr arg)))))
@

The function \verb+indexed+ produces a one-dimensional list of \verb+gname+s.
The inputs are the string \verb+anamestring+, and the numbers \verb+starting+
and \verb+ending+.  For convenience, \verb+anamestring+ may be a \verb+gname+,
in which case, the first component of its name is used.  A precondition is that
\verb+starting+ be less than or equal to \verb+ending+. 
Here is the function itself.  
<<>>=
(defun indexed (anamestring starting ending)
  (if (consp anamestring)
      (indexed (car anamestring) starting ending)
    (mapcar #'(lambda (x) 
               (@ anamestring x))
            (seq starting ending))))
@

The function
\verb+vecname+ returns \verb+NIL+ if the argument
\verb+aname+ has no brackets in it, and otherwise it returns 
the stuff in front of the bracket.  
<<>>=
(defun vecname (aname)
  (let* ((thename (cond ((consp aname) (car aname))
                        ((symbolp aname) (string aname))
                        (t aname)))
         (loc (search "[" (string thename))))
    (cond ((null loc) NIL)
           (t (subseq (string thename) 0 loc)))))
(defun vec-index (aname)
  (let* ((thename (cond ((consp aname) (car aname))
                        ((symbolp aname) (string aname))
                        (t aname)))
         (ss (string thename))
         (loc (search "[" ss))
         (loc2 (search "]" ss)))
    (cond ((null loc) NIL)
           (t (read-from-string (subseq (string thename) (+ 1 loc) loc2))))))
@

\subsection*{Data Structures}
A \verb+definition+ is something that requires an
S-expression to be computed.
The field \verb+sform+ should contain a single LISP S-expression.
The field \verb+tag-flag+ is true if the \verb+sform+ is a list,
sum, or product of atomic items, and false otherwise.
<<>>=
(defstruct definition
  name
  tag-flag
  involves-transitions
  sform)
(defun new-definition (&key name sform (tag NIL) (involves-transitions NIL))
  (make-definition :name name 
                   :sform sform 
                   :involves-transitions involves-transitions
                   :tag-flag tag))
(defun copy-definition (adef)
  (make-definition :name (definition-name adef)
                   :sform (copy-tree (definition-sform adef))
                   :tag-flag (definition-tag-flag adef)))
@

The function \verb+def-name-in+ returns \verb+T+ if the name of the
argument \verb+d1+ is in the name of \verb+d2+.  Both \verb+d1+ and
\verb+d2+ must be definitions.  A simple example:
\begin{verbatim}
> (setf z1 (new-definition :name $s1 :sform 2))
#S(DEFINITION NAME #S(GNAME NAME (S1)) TAG-FLAG NIL SFORM 2)
> (setf z2 (new-definition :name $(s1 k1) :sform 3))
#S(DEFINITION NAME #S(GNAME NAME (K1 S1)) TAG-FLAG NIL SFORM 3)
> (def-name-in z1 z3)
T
> (def-name-in z3 z1)
NIL
\end{verbatim}
Here is the function itself:
<<>>=
(defun def-name-in (d1 d2)
  (name-included-p (definition-name d1) (definition-name d2)))
@

The function \verb+same-def-p+ returns \verb+T+ if the definitions have
the same name, and \verb+NIL+ otherwise.
<<>>=
(defun same-def-p (d1 d2)
  (eq (definition-name d1) (definition-name d2)))  ; same-name-p to eq
@

A \verb+transition+ object contains a \verb+name+, a symbolic
expression \verb+value+, and a 
\verb+from+ and \verb+to+ state.
<<>>=
(defstruct transition
  name
  value
  from
  to
  tags)
@

<<>>=
(defun gen-transition-name (from to)
  (new-gname
   (intern (concatenate 'string (to-string (gensym))
                                "_from" (to-string from) 
                                "_to" (to-string to)))))
(defun new-transition (&key from to value (tag NIL))
  (make-transition :name (gen-transition-name from to)
                   :from from
                   :to to
                   :value value
                   :tags tag))
(defun deep-copy-transition (tr1)
  (make-transition :name (transition-name tr1)
                   :from (transition-from tr1)
                   :to (transition-to tr1)
                   :value (copy-tree (transition-value tr1))
                   :tags (copy-tree (transition-tags tr1))))
@

We also have four more utility
functions.  The first argument to \verb+is-from?+ is a
transition; the second is a state; the function returns \verb+T+ if
the transition is \verb+from+ the state.  The function \verb+is-to?+ is
analogous.  We also have \verb+is-from-any?+ which accepts a transition
and a list of states; if the transition is from any of the states in 
the list, the function returns non-\verb+NIL+, and analogously for 
\verb+is-to-any?+.  
<<>>=
(defun is-from (t1 s1)
  (eq (transition-from t1) s1))   ; eq was same-name-p
(defun is-to (t1 s1)
  (eq (transition-to t1) s1))     ; eq was same-name-p
(defun is-from-any (t1 slist)
  (member-if #'(lambda (x) (is-from t1 x)) slist)) 
(defun is-to-any (t1 slist)
  (member-if #'(lambda (x) (is-to t1 x)) slist)) 
@

First, we define the \verb+model+ data structure.  The fields \verb+states+,
\verb+history+,
\verb+parameters+, \verb+transitions+, 
\verb+innerfunctions+ and \verb+definitions+ 
each contain a list of the appropriate type of objects.
<<>>=
(defstruct model
  name
  history
  states
  parameters
  transitions
  definitions
  innerfunctions)
(defun new-model (name)
  (make-model :name (new-gname name)
              :history ()
              :states (make-hash-table :test #'equal)
              :parameters (make-hash-table :test #'equal)
              :transitions ()
              :innerfunctions (make-hash-table :test #'equal)
              :definitions ()))
(defun copy-hash-table (ht)
  (let ((newtbl (make-hash-table :test (hash-table-test ht))))
    (maphash #'(lambda (key value)
                       (setf (gethash key newtbl) value))
             ht)
    newtbl))
(defun deep-copy-model (amodel)
  (make-model :name (model-name amodel)
              :history ()
              :states (copy-hash-table (model-states amodel))
              :parameters (copy-hash-table (model-parameters amodel))
              :transitions (mapcar #'deep-copy-transition (model-transitions amodel))
              :innerfunctions (copy-hash-table (model-innerfunctions amodel))
              :definitions (mapcar #'deep-copy-definition (model-definitions amodel))))
@

The function \verb+defs-matching1+ returns all definitions in the
argument \verb+amodel+ that include the name of the definition \verb+adef+.
The function \verb+defs-matching+ returns all definitions in the argument
\verb+amodel+ that match anything in \verb+adeflist+.  
<<>>=
(defun defs-matching1 (amodel adef)
  (remove-if-not #'(lambda (x) (def-name-in adef x)) 
                 (model-definitions amodel)))
(defun defs-matching (amodel adeflist)
  (remove-duplicates (defs-matching-aux amodel adeflist NIL)
                     :test #'same-def-p))
(defun defs-matching-aux (amodel adeflist sofar)
  (cond ((null adeflist) sofar)
        (t (defs-matching-aux amodel 
                              (cdr adeflist)
                              (append (defs-matching1 model (car adeflist))
                                      sofar)))))
@

First, the use of \verb+defstruct+ has already generated 
\verb+model-states+, \verb+model-parameters+, \verb+model-name+, 
\verb+model-transitions+, and \verb+model-definitions+.  

We have the following \verb+add+ functions, which append an object
of the indicated type to the appropriate field.
<<>>=
(defun add-state-nocheck (amodel thestate)
  (setf (gethash thestate (model-states amodel)) (symbolize thestate)))
(defun add-state (amodel thestate)
  (add-state-nocheck amodel thestate))
(defun add-state-list (amodel alist)
  (dolist (nn alist)
    (add-state amodel nn)))
(defun remove-state (amodel astate)
  (remhash astate (model-states amodel))
  (setf (model-transitions amodel)    ; eq was same-name-p
        (remove-if #'(lambda (x) (or (eq (transition-from x) astate)
                                     (eq (transition-to x) astate)))
                   (model-transitions amodel)))
  (let* ((defs (model-definitions amodel))
         (hdr NIL))
    (dolist (dd defs)
      (setf hdr (car (definition-sform dd)))
      (setf (definition-sform dd)              ; eq was same-name-p
            (cons hdr (remove-if #'(lambda (x) (eq x astate)) 
                                 (cdr (definition-sform dd))))))))
@

<<>>=
(defun add-parameter-nocheck (amodel theparameter)
  (setf (gethash theparameter (model-parameters amodel)) (symbolize theparameter)))
(defun add-parameter (amodel theparameter)
  (add-parameter-nocheck amodel theparameter))
(defun add-parameter-list (amodel alist)
  (dolist (nn alist)
    (add-parameter amodel nn)))
(defun remove-parameter (amodel aparam)
  (remhash aparam (model-parameters amodel)))
@

<<>>=
(defun add-definition (amodel aname &key value (tag NIL) (involves-transitions NIL) (add-tag NIL))
  (let* ((vecpart (vecname aname))
	 (addtag (if (not add-tag) NIL (not (null vecpart))))
         (thename (if (gname-p aname) aname (new-gname aname)))
         (thedefinition (new-definition :name thename
                                        :sform value
                                        :involves-transitions involves-transitions
                                        :tag tag)))
    (if addtag
        (set-tag amodel :targets (list (definition-name thedefinition)) 
                        :tag (new-gname vecpart)))
                                      ; eq was same-name-p
    (let* ((defs (model-definitions amodel)))
      (setf (model-definitions amodel) 
            (cons thedefinition defs)))))
(defun insert-definition (amodel adef)
  (setf (model-definitions amodel) (cons adef (model-definitions amodel))))
@

<<>>=
(defun redefine (amodel aname &key value (add-tag T))
  (let* ((thename (new-gname aname)))
    (remhash aname (model-parameters amodel))
    (add-definition amodel aname :value value :add-tag add-tag)))
(defun remove-definition (amodel aname)
  (let ((thename (new-gname aname)))
    (setf (model-definitions amodel)
          (remove-if #'(lambda (x) (same-name-p thename (definition-name x)))
                     (model-definitions amodel)))))
(defun clear-all-definitions (amodel)
  (let ((thedefs (model-definitions amodel)))
    (setf (model-definitions amodel) NIL)
    thedefs))
@

<<>>=
(defun add-transition (amodel &key from to rate (tag NIL))
  (let ((thetransition (new-transition :from from
                                       :to to
                                       :value rate
                                       :tag tag)))
    (setf (model-transitions amodel) 
          (cons thetransition (model-transitions amodel)))
    thetransition))
(defun remove-transition (amodel aname)
  (let ((thename (new-gname aname)))
    (setf (model-transitions amodel)
          (remove-if #'(lambda (x) (same-name-p thename (transition-name x)))
                     (model-transitions amodel)))))
(defun change-transition-to (amodel aname new-to-state)
  (let* ((thename (new-gname aname))
         (thetr (find-if #'(lambda (x) (same-name-p thename (transition-name x)))
                        (model-transitions amodel))))
    (setf (transition-to thetr) (new-gname new-to-state))))
(defun clear-all-transitions (amodel)
  (let ((alltrs (model-transitions amodel)))
    (setf (model-transitions amodel) NIL)
    alltrs))
@

<<>>=
(defun transitions-from (mod st)
  (remove-if-not #'(lambda (tt) (is-from tt st)) (model-transitions mod)))
(defun transitions-to (mod st)
  (remove-if-not #'(lambda (tt) (is-to tt st)) (model-transitions mod)))
(defun transitions-from-any (mod st)
  (remove-if-not #'(lambda (tt) (is-from-any tt st)) (model-transitions mod)))
(defun transitions-to-any (mod st)
  (remove-if-not #'(lambda (tt) (is-to-any tt st)) (model-transitions mod)))
@

Now, tags, implemented as definitions.
<<>>=
(defun get-dname (x)
  (cond ((definition-p x) (definition-name x))
        (t x)))
(defun set-tag (amodel &key targets tag)
  (if (null tag) (error "set-tag: null tag"))
  (let* ((targets1 (if (not (gname-p targets)) targets (list targets)))
         (tagname (get-dname tag))
         (defnames (mapcar #'definition-name (model-definitions amodel)))
         (tagloc (position-if        ; same-name-p changed to eq
                      #'(lambda (x) (eq tagname x))
                      defnames))
         (targlist (cond ((consp targets1) targets1)
                         (t (list targets1)))) ; huh?
         (trs (if (transition-p (car targlist)) T NIL)))
    (cond ((null tagloc)
           (add-definition amodel tag :value (cons 'LIST targlist) :tag T :involves-transitions trs))
          (t
           (let* ((curtag (nth tagloc (model-definitions amodel)))
                  (ds (definition-sform curtag))
                  (otargs 
                   (remove-duplicates
                    (append targlist (copy-list (cdr ds)))
                    :test #'same-object-name-p)))
             (setf (definition-sform curtag)
                   (cons 'LIST otargs)))))
    (dolist (tg targlist)
      (if (transition-p tg) 
          (setf (transition-tags tg) (cons tag (transition-tags tg)))))))
(defun show-tag (amodel atag)
  (let* ((thename (new-gname atag))   ; same-name-p made to eq
         (dm (find-if #'(lambda (x) (eq thename (definition-name x)))
                      (model-definitions amodel))))
    (if (not (null dm))
        (definition-sform dm))))
@

Here is a shortcut constructor.
<<>>=
(defun add-all-flows (amodel &key states parameter)
  (add-state-list amodel states)
  (let ((pname-prefix (string (show parameter)))
        (len (length states)))
    (dotimes (ii len)
      (dotimes (jj len)
        (if (not (equal ii jj))
          (let ((curpar (@ pname-prefix (list jj ii)))
                (fromstate (nth ii states))
                (tostate (nth jj states)))
            (add-parameter amodel (new-gname curpar))
            (add-transition amodel :from fromstate :to tostate
                            :rate (list '* fromstate curpar))))))))
@

<<>>=
; should be more general--use an expression not just a parameter
(defun add-flow-series (amodel &key states parameter (tag NIL) (stratify T))
  (add-state-list amodel states)
  (if tag
      (set-tag amodel :targets states :tag tag))
  (if (null stratify)
      (progn
       (add-parameter amodel parameter)
       (let ((len (length states)))
         (dotimes (ii (- len 1))
           (let ((fromstate (nth ii states))
                 (tostate (nth (+ 1 ii) states)))
             (add-transition amodel :from fromstate :to tostate
                             :rate (list '* fromstate parameter))))))
    (progn
      (let ((pname-prefix (string (show parameter)))
            (len (length states)))
        (dotimes (ii (- len 1))
          (let ((curpar (@ pname-prefix ii))
                (fromstate (nth ii states))
                (tostate (nth (+ 1 ii) states)))
            (add-parameter amodel (new-gname curpar))
            (add-transition amodel :from fromstate :to tostate
                            :rate (list '* fromstate curpar))))))))
@

An important pattern consists of having a stratified parameter 
associated with every state in a model.  There are several such
patterns: the first is to have a flow from each state into another
state (for example, stage specific disease mortality).  Another such
pattern is to have a cost rate associated with each state, and to accrue total
costs according to these costs.  The first of these patterns is handled
by \verb+add-all-flows-into+; this might be a job for \verb+add-transition+
in the future.
<<>>=
(defun add-all-flows-into (amodel &key from to parameter (stratify T))
  (let ((np
         (if stratify
             (mapcar #'(lambda (x) (merge-names x parameter)) from)
           (dup parameter (length from)))))
    (if stratify
        (dolist (pp np)
          (add-parameter amodel pp))
      (add-parameter amodel parameter))
    (dotimes (ii (length from))
      (add-transition amodel :from (nth ii from)
                             :to to
                             :rate (list '* (nth ii from) (nth ii np)))))) 
@

The second such pattern is one which will simply generate the list of
stratified parameter names.  The list of parameter names is returned.
<<>>=
(defun add-parameters-by-states (amodel states param)
  (let ((ps (mapcar #'(lambda (x) (merge-names x param)) states)))
    (dolist (pp ps)
      (add-parameter amodel pp))
    ps))
@

<<>>=
(defun dropvec1 (astring)
  (let* ((aa (vecname astring)))
    (if (null aa)
        (if (symbolp astring) astring (intern astring))
      (intern aa))))
(defun dropvec (aname)
  (new-gname (mapcar #'dropvec1 aname)))
(defun name-contains (a1 a2)
  (name-included-p a1 (dropvec (get a2 'namelist))))
(defun get-parameter (amodel aname)   ; same-name-p changed to eq
  (remove-if-not #'(lambda (x) (eq aname x))
    (model-parameters amodel)))
(defun gpc-aux (curlist arg)
  (cond ((null arg) curlist)
        ((null curlist) NIL)
        (t (gpc-aux
            (remove-if-not #'(lambda (x) (name-contains (car arg) x)) 
                           curlist)
            (cdr arg)))))
; this is modified from the CLTL web site.
(defun collect (bool &optional key value)
  (list bool key value))
(defun ht-filter (ht pred)
  (with-hash-table-iterator (getnext ht) 
    (do ((cur (multiple-value-call #'collect (getnext))  
              (multiple-value-call #'collect (getnext)))
         (col () (if (and (not (null (car cur))) 
                          (funcall pred (cadr cur)))
                       (cons (cadr cur) col)
                     col)))
      ((null (car cur)) col))))
(defun states-without (amodel aname)    ; same-name-p => eq
  (ht-filter (model-states amodel) #'(lambda (x) (eq aname x))))
(defun parameters-without (amodel aname)    ; same-name-p => eq
  (ht-filter (model-parameters amodel) #'(lambda (x) (eq aname x))))
(defun st (amodel)
  (ht-filter (model-states amodel) #'(lambda (x) T)))
(defun pa (amodel)
  (ht-filter (model-parameters amodel) #'(lambda (x) T)))
@

This defines \verb+name+ to be the sum of the states times the 
parameters.
<<>>=
(defun add-definition-linsum (amodel &key name states parameters)
  (let ((ex (cons '+ (mapcar #'(lambda (x y) (list '* x y)) states parameters))))
    (add-definition amodel name :value ex)))
@

\subsection*{Cross-stratification}

<<>>=
(defun flatten-one-level (alist)
  (labels ((add-to-ans (ulist vlist)
             (do* ((ans0 vlist (cons (car cur0) ans0))
                   (cur0 ulist (cdr cur0)))
                 ((null cur0) ans0))))
    (do* ((ans () (add-to-ans (car cur) ans))
          (cur alist (cdr cur)))
        ((null cur) (reverse ans)))))
@

<<>>=
(defun outermap (fn alist blist)
  (mapcar 
   #'(lambda (y) (mapcar #'(lambda (x) (funcall fn x y)) alist))
   blist))
@

See Graham, page 42.
<<>>=
; treacherous: any symbol with a property list gets counted as 
;   a gname
(defun gname-p (xx)
  (and (symbolp xx) (get xx 'namelist)))
(defun stratify-name (pp ss nd production-function)
  (let ((the-assoc (gethash pp nd)))
    (cond ((null the-assoc)
           (funcall production-function pp ss))
          ((eq the-assoc 0) pp)
          (t 
           (funcall the-assoc pp ss)))))
(defun stratify-expression (newname tree nodep production-function)
  (if (gname-p tree)    
      (stratify-name tree newname nodep production-function)
    (if (atom tree)
        tree
      (cons (stratify-expression newname (car tree) nodep production-function)
            (stratify-expression newname (cdr tree) nodep production-function)))))
@

Note the exclusion criteria for transitions. Self-transitions sometimes
arise in crossing one model by another when you use a 
``production function'' (combining rule for names) that is something
other than merge-names, the default.  Also excluded are transitions that connect
two states that were both generated by nondefault combining rules
(production functions other than merge-names); this keeps you from
accidentally duplicating transitions in the part of the model you did
not stratify.
<<>>=
(defun stratify-transition (tr1 s1 no-dependence production-function)
  (multiple-value-bind (newfrom from-p) 
                       (funcall production-function (transition-from tr1) s1)
    (multiple-value-bind (newto to-p) 
                         (funcall production-function (transition-to tr1) s1)
      (let ((newvalue (stratify-expression s1 
                                           (transition-value tr1)
                                           no-dependence
                                           production-function)))
        (if (or (equal newfrom newto)
                (and (not from-p) (not to-p))
                (and (source-p newfrom) (not to-p))
                (and (not from-p) (sink-p newto)))
            NIL
          (new-transition :from newfrom
                          :to newto
                          :value newvalue
                          :tag (transition-tags tr1)))))))
@

<<>>=
(let ((ctr 0))
  (defun gen-tmp-dname ()
    (let ((symst (concatenate 'string "_DEF" (write-to-string ctr))))
      (incf ctr)
      (new-gname symst))))
(defun pdm (x m)
  (position (definition-name x) m :test #'equalp))
(defun add-lin-trans (amodel mset pname)
  (let* ((dsp (split-list 
               #'(lambda (x) (any-included-p mset (definition-name x))) 
               (model-definitions amodel)))
         (dy0 (copy-list (car dsp)))
         (dy (sort dy0 #'(lambda (x y) (<= (pdm x m)
                                           (pdm y m)))))
         (ny (cadr dsp))
         (cdy (mapcar 
               #'(lambda (x) (new-definition :name (definition-name x)
                                             :sform NIL))
               dy))
         (dy1 (progn
                (dolist (dd dy)
                  (setf (definition-name dd) (gen-tmp-dname)))
                dy))
         (newdnames (mapcar #'definition-name dy1))
         (ns (seq 0 (- (length newdnames) 1))))
    (dotimes (ii (length cdy))
      (setf (definition-sform (nth ii cdy))
        (cons '+
              (mapcar #'(lambda (x y) (list '* x y))
                      (mapcar #'(lambda (x) (@ pname (list ii x))) ns)
                      newdnames))))
    (setf (model-definitions amodel)
          (sort (copy-list (append (model-definitions amodel) cdy)) #'def<))
    (mapcar #'definition-name cdy)))
@

<<>>=
(defun select (aa bool)
  (cond ((equal (length aa) (length bool))
         (let ((ans NIL))
           (dotimes (ii (length aa))
             (if (not (null (nth ii bool)))
                 (setf ans (cons (nth ii aa) ans))))
           ans))
        (t (error "lists not equal"))))
@

Reverse all boolean items in a boolean list:
<<>>=
(defun negate (bool)
  (let ((ans ()))
    (dolist (bb bool)
      (setf ans (cons (not bb) ans)))
    (reverse ans)))
@

<<>>=
(defun as-matrix (&key definition matrix-name)
  (cond ((and (definition-p definition) (gname-p matrix-name))
         (list definition matrix-name))
        (t (error "as-matrix: wrong types."))))
@

We need to be able to remove duplicates without testing all
elements pairwise.
<<>>=
(defun l< (aa bb)
  (string< (symbolize aa) (symbolize bb)))
(defun gs-remove-duplicates (xx fn)
  (let* ((thelist (sort (copy-list xx) #'l<)))
    (do ((cur thelist (cdr cur))
         (ans () (if (funcall fn (car ans) (car cur))
                     ans
                   (cons (car cur) ans))))
      ((null cur) ans))))
@

The function \verb+stratify-group+ stratifies each member of 
\verb+pars+ by every member of \verb+states+, excluding \verb+nodep+,
and returns a list.
<<>>=
(defun stratify-group (pars nodep states production-function par-ht)
  (let* ((the-assoc ())
         (par-stratified ()))
    (dolist (pp pars)
      (setf the-assoc (gethash pp nodep))  
      (dolist (ss states)    ; the eq below was a same-name-p
        (setf (gethash (stratify-name pp ss nodep production-function) par-ht) T))))
  par-ht)
@
It is worth remembering that in some implementations of Lisp, the 
\verb+remove-duplicates+ function compares things pairwise instead of sorting.

Sometimes you want to take definitions and stratify only the
states in them.  For instance, take a definition whose expression
is \verb+(LIST X Y)+ and stratify by \verb+A1+ and \verb+A2+, yielding
\verb+(LIST A1-X A1-Y A2-X A2-Y)+, rather than two definitions
with \verb+(LIST A1-X A2-X)+ and \verb+(LIST A1-Y A2-Y)+.  
<<>>=
(defun stratify-defs (defs no-dependence states production-function)
  (let* ((def-stratified ())
         (the-assoc ())
         (doneit NIL))
    (dolist (pp defs)
      (setf doneit NIL)
      (setf the-assoc (gethash (definition-name pp) no-dependence))
      (dolist (ss states)   ; the eq was a same-name-p
        (cond ((null the-assoc)
               (setf def-stratified
                    (cons (stratify-definition pp ss no-dependence production-function) def-stratified)))
              (t (if (not doneit)
                     (progn
                      (setf def-stratified (cons pp def-stratified))
                      (setf doneit T)))))))
    def-stratified))
(defun stratify-definition (def1 s1 no-dependence production-function)
  (if (definition-involves-transitions def1)
    NIL
    (let* ((dsf (definition-sform def1))
           (newvalue 
            (if (and (consp dsf) (equal (car dsf) 'LIST))
                NIL
              (stratify-expression s1 dsf no-dependence production-function))))
      (new-definition :name (merge-names (definition-name def1) s1)
                      :sform newvalue
                      :tag (definition-tag-flag def1)))))
(defun regenerate-tags (thedefs)
  (let ((taghash (make-hash-table :test #'equal)))
    (dolist (dd thedefs)
      (if (definition-tag-flag dd)
          (setf (gethash (definition-tag-flag dd) taghash) 
                (cons (definition-sform dd)
                      (gethash (definition-tag-flag dd) taghash)))))
    (maphash #'(lambda (k v)
                (new-definition :name k
                                :sform (cons '+ v)
                                :tag NIL))
             taghash)))
@

In stratifying, we must apply a function to the parameter or state
(say A)
of the first model, and the state B of the second, to produce the
new parameter or state.  We use the following conventions.  The
necessary information is stored in a hash table under the name
of the parameter or state of the first model.  
If nothing is stored under A, then it is assumed that when combined
with B, the result is simply AB, i.e. complete stratification. If
the number 0 is stored, we assume that the result is simply A, no
stratification at all.  Finally, if a lambda expression is stored
under A, this expression must return the name of the state or
parameter which results.  This mechanism lets us choose fewer
new parameters than a full stratification would lead to.
<<>>=
(defun new-dictionary ()
  (make-hash-table :test #'equal))
(defun set-dependence (ht &key of to)
  (cond ((consp of)
         (dolist (pp of)
           (setf (gethash pp ht) to)))
        (t (setf (gethash of ht) to)))
  ht)
(defun suppress-dependence (ht &key of)
  (set-dependence ht :of of :to 0))
(defun set-combining-function (ht &key target level)
  (set-dependence ht :of target :to level))
@

Next, handle tags made from transitions.  Collect up all the transitions
with given tags.  Why is this done here? So it won't have to be done at
compile time.
<<>>=
(defun gen-transtags (translist)
  (let ((tbl (make-hash-table :test #'equal)))
    (dolist (tr translist)
      (if (not (null (transition-tags tr)))
          (dolist (tt (transition-tags tr))
            (if (null (gethash tt tbl))
                (setf (gethash tt tbl) (list (transition-name tr)))
              (setf (gethash tt tbl) (cons (transition-name tr) (gethash tt tbl)))))))
    (let ((ans NIL))
      (maphash #'(lambda (x y) (setf ans (cons (new-definition :name x :sform (cons 'LIST y) :tag T :involves-transitions T) ans))) tbl)
      ans)))
@

Next, the cross-stratification function.  
There are several issues.  First, you may specify that certain parameters
and definitions may not be stratified by certain states; this is
done using \verb+no-dependence+ which is a hash table.
<<>>=
(defun hash-it-all (alist atest)
  (let ((newht (make-hash-table :test atest)))
    (dolist (cur alist)
      (setf (gethash cur newht) T))
    newht))
; taken from pleac.sourceforge.net on 7 jan 2015
(defun merge-hash-tables (&rest tables)
  (let ((union
         (make-hash-table
          :test (first
                 (sort (mapcar #'hash-table-test tables) #'>
                       :key (lambda (test)
                              (ecase test
                                (eq 0)
                                (eql 1)
                                (equal 2)
                                (equalp 3)))))
          :size (reduce #'max (mapcar #'hash-table-size tables)))))
    (dolist (table tables)
      (maphash (lambda (key val) (setf (gethash key union) val)) table))
    union))
; to get usual cartesian behavior, make
;   production function equal to #'merge-names-fn
; the production-function should return two values; first the new name,
; second, t if the name is new, NIL else
(defun cross (model-a model-b xclude-dict &optional (production-function #'merge-names-fn)  (verbose NIL))
  (let* ((nd (if (null xclude-dict) 
                 (make-hash-table)
               xclude-dict)) 
         ; first, start the new model
         (modelnames (merge-names (model-name model-a)
                                  (model-name model-b)))
         (newpar-ht (make-hash-table :test #'equal))
         (zzz (if verbose 
                  (format t "First call to stratify-group~%")
                NIL)) 
         ; step 1: take care of the parameters
         (ignore1 (stratify-group (pa model-a)
                                  nd
                                  (st model-b)
                                  production-function
                                  newpar-ht))
         (zzz1 (if verbose 
                  (format t "Second call to stratify-group~%")
                NIL)) 
         (ignore2 (stratify-group (pa model-b)
                                  nd
                                  (st model-a)
                                  production-function
                                  newpar-ht))
         ; step 2: take care of the definitions. 
         (zzz2 (if verbose 
                  (format t "First call to stratify-defs~%")
                NIL)) 
         (def-a-stratified1 (stratify-defs (model-definitions model-a)
                                           nd
                                           (st model-b)
                                           production-function))
         (zzz3 (if verbose 
                  (format t "Second call to stratify-defs~%")
                NIL)) 
         (def-b-stratified1 (stratify-defs (model-definitions model-b)
                                           nd
                                           (st model-a)
                                           production-function))
         (zzka (if verbose
                   (format t "done with stratify-defs~%")
                NIL))
         (def-a-stratified (remove-if #'null def-a-stratified1))
         (zzka1 (if verbose
                   (format t "after def-a-stratified~%")
                NIL))
         (def-b-stratified (remove-if #'null def-b-stratified1))
         (zzka2 (if verbose
                   (format t "after def-b-stratified~%")
                NIL))
         ; step 3: take care of tag definitions
         (hdefs-a (regenerate-tags (append def-a-stratified def-b-stratified)))
         (zzka3 (if verbose
                   (format t "after hdefs-a")
                NIL))
         ; step 4: handle transitions
         (zzz4 (if verbose 
                  (format t "first call to stratify-transitions~%")
                NIL)) 
         (t1-stratified
          (remove-if #'null 
           (flatten-one-level
            (outermap #'(lambda (d s) 
                         (stratify-transition d s nd production-function)) 
                      (model-transitions model-a)
                      (st model-b)))))
         (zzz5 (if verbose 
                  (format t "second call to stratify-transitions~%")
                NIL)) 
         (t2-stratified
          (remove-if #'null 
           (flatten-one-level
            (outermap #'(lambda (d s) 
                         (stratify-transition d s nd production-function))
                      (model-transitions model-b)
                      (st model-a)))))
         (modeltransitions (append t1-stratified t2-stratified))
         ; handle definitions involving transitions
         (transtags (gen-transtags modeltransitions))
         (zzz6 (if verbose 
                  (format t "outer map~%")
                NIL)) 
         (modelstates 
          (flatten-one-level 
           (outermap #'(lambda (x y) (funcall production-function x y))
                     (st model-a)
                     (st model-b))))
         (zzz7 (if verbose 
                  (format t "definitions~%")
                NIL)) 
         (modeldefinitions 
          (append def-a-stratified def-b-stratified hdefs-a transtags))
         (modelinnerfunctions
          (merge-hash-tables (model-innerfunctions model-a) (model-innerfunctions model-b))))
  (make-model :name modelnames
              :history ()
              :states (hash-it-all modelstates #'equal)
              :parameters newpar-ht
              :transitions modeltransitions
              :innerfunctions modelinnerfunctions
              :definitions modeldefinitions)))
@

\section*{Compilation to Ordinary Differential Equations}
We need to perform a topological sort on the list of definitions and
transitions.  First, some functions that maintain definitions and
transitions in a list; the definition or transition itself is \verb+item+
and the dependencies are listed as \verb+dependency-list+ at first.  Only
other definitions or transitions are listed there.  This is used as a
temporary variable for the sorting.

<<>>=
(defstruct obj
  item
  dependency-list)
(defun flatten (tree)
  (cond ((null tree) NIL)
        ((atom tree) (list tree))
        (t (flatten-aux tree '()))))
; check this function! note the orders get scrambled a bit
(defun flatten-aux (curtree curans)
  (let* ((newtree curtree)
         (newans ()))
    (tagbody
     the-top
     (cond ((null newtree) (go the-end))
           ((not (consp (car newtree)))
            (progn
             (setf newans (cons (car newtree) newans))
             (setf newtree (cdr newtree))
             (go the-top)))
           ((null (cdr newtree))
            (setf newans (flatten-aux (car newtree) newans)))
           (t (setf newans (append (flatten-aux (car newtree) '())
                                   (flatten-aux (cdr newtree) newans)))))
     the-end)
     (append (reverse newans) curans)))
; the next works well on small functions but overflows the stack on big ones 
(defvar *maxapply* 60)
(defun fappend (flist)
  (if (< (length flist) 60)
      (apply #'append flist)
    (fappend-aux flist (length flist) '())))
(defun fappend-aux (flist nn sofar)
  (if (< nn 60)
      (append (apply #'append flist) sofar)
    (fappend-aux (cdr flist) (- nn 1) (append (car flist) sofar))))
(defun getname (u)
  (cond ((gname-p u) u)         
        ((definition-p u) u)
        ((transition-p u) u)
        (t NIL)))

(defun allnames (xp)
  (ht-filter (allnames-aux (oform xp) (make-hash-table :test #'equal)) #'(lambda (x) T)))
(defun allnames-aux (cur ht)
  (cond ((null cur) ht)
        ((not (listp cur))
         (progn
          (setf ztemp (getname cur))
          (if (not (null ztemp))
            (setf (gethash ztemp ht) T)))
         ht)
        (t (dolist (c1 cur)
             (setf ht (allnames-aux c1 ht)))
           ht)))

(defun same-object-name-p (x y)    ; the eq was a same-name-p
  (eq (getname x) (getname y)))
(defun oname (o1)
  (cond ((gname-p (obj-item o1)) (obj-item o1))   
        ((definition-p (obj-item o1)) (definition-name (obj-item o1)))
        ((transition-p (obj-item o1)) (transition-name (obj-item o1)))
        (t NIL)))
(defun oform (o1)
  (cond ((gname-p o1) NIL)   
        ((definition-p o1) (definition-sform o1))
        ((transition-p o1) (transition-value o1))
        (t NIL)))
(defun gen-obj (anobj)
  (make-obj :item anobj 
            :dependency-list (allnames anobj)))
(defun nodep (aobj)
  (null (obj-dependency-list aobj)))
(defun def-is-list-p (x)
  (if (consp (definition-sform x))
      (equalp 'list (car (definition-sform x)))
    NIL))
@

Don't change gen-dlist to anything more "elegant". You get
invocation stack history overflows on huge models if you
aren't careful.  Leave the tagbody's. It's fine.
<<>>=
(defun gen-dlist (amodel)
  (let ((ans ())
        (g 0)
        (cur (st amodel)))
    (tagbody
     the-top
     ;(prin1 (incf g))
     (setf ans (cons (gen-obj (car cur)) ans))
     (setf cur (cdr cur))
     (if (not (null cur))
         (go the-top)))  
    (setf cur (pa amodel))
    (prin1 "parameters...")
    (tagbody
     the-top
     ;(prin1 (incf g))
     (setf ans (cons (gen-obj (car cur)) ans))
     (setf cur (cdr cur))
     (if (not (null cur))
         (go the-top)))
    (setf cur (model-transitions amodel))
    (prin1 "transitions...")
    (tagbody
     the-top
     ;(prin1 (incf g))
     (setf ans (cons (gen-obj (car cur)) ans))
     (setf cur (cdr cur))
     (if (not (null cur))
         (go the-top)))
    (setf cur (remove-if #'def-is-list-p
                         (model-definitions amodel)))
    (prin1 "definitions...")
    (prin1 (length cur))
    (terpri)
    (tagbody
     the-top
     ;(prin1 (incf g))
     ;(terpri)
     ;(system "date")
     ;(terpri)
     (setf ans (cons (gen-obj (car cur)) ans))
     (setf cur (cdr cur))
     (if (not (null cur))
         (go the-top)))
     ans)) 
@

<<>>=
(defun remove-dependencies (obj htbl)
  (setf (obj-dependency-list obj)
        (remove-all htbl
                    (obj-dependency-list obj) 
                    '())))
@

<<>>=
(defun cie (x) 
  (if (consp x) (car x) NIL))
(defun rsp (x)
  (remove (coerce " " 'character) x))
(defun remove-all (htbl curseq sofar)
  (let* ((newseq curseq)
         (newsofar sofar))
    (tagbody
     the-top
     (cond ((null newseq) (go the-end))
           ((and (not (null (car newseq))) (gethash (show (cie newseq)) htbl))
            (setf newseq (cdr newseq))
            (go the-top))
           (t 
            (setf newsofar (cons (car newseq) newsofar))
            (setf newseq (cdr newseq))
            (go the-top)))
     the-end)
     newsofar))
@

<<>>=
; all objects in deps-to-rem must be removed from the dependency
;   list of each object in the list  objs
(defun remove-all-dependencies (objs deps-to-rem)
  (let* ((ht (make-hash-table :test #'equal)))
    (dolist (dd deps-to-rem)
      (setf (gethash (show (oname dd)) ht) T))
    (dolist (oo objs)
      (remove-dependencies oo ht))
    objs))
@

<<>>=
(defun dep-split-aux (uns srt)
  (let* ((unsorted uns)
         (sorted srt)
         (curlen 0)
         (prevlen 0)
         (first T)
         (s1 ())
         (new-nodep ())
         (new-dep ()))
    (tagbody
      the-top
      (if first
          (progn 
           (setf prevlen (length unsorted))
           (setf first NIL))
        (progn
         (if (= curlen prevlen)
             (progn 
              (setf *sorted* sorted) (setf *unsorted* unsorted)
              (error "*** Compile error occurred due to a cycle or missing declaration."))
           (progn
            (setf prevlen curlen)
            (setf curlen (length unsorted))))))
      (format t "compiling: (length unsorted) ~A" (length unsorted))
      (terpri)
      (setf s1 (split-list #'nodep unsorted))
      (setf sorted (append sorted (car s1)))
      (setf new-dep (cadr s1))
      (setf unsorted (remove-all-dependencies new-dep sorted))
      (if (not (null unsorted))
          (go the-top)))
    sorted))
@

<<>>=
(defun topo-sort (amodel)
  (let* ((amodel2 (copy-tree amodel))
         (g1 (gen-dlist amodel2))
         (zzzk (prin1 "done with dependency list")))
    (mapcar #'obj-item
            (dep-split-aux g1 NIL))))
@

<<>>=
(defun dtfilter (oo)
  (cond ((definition-p oo) oo)
        ((transition-p oo) oo)
        (t NIL)))
(defun genname (i s)
  (intern (concatenate 'string (string s) (format NIL "~A" i))))
(defun assym1 (a)
  (intern (remove " " (show a) :test #'string=)))
(defun to-def (asym parlist)
  (let ((outs ()))
    (dotimes (ii (length parlist))
      (setf outs (cons (list (assym1 (nth ii parlist))
                             (list 'AREF asym ii))  
                       outs)))
    outs))
(defun gen-index (statelist)
  (let* ((ht (make-hash-table :test #'equalp)))
    (dotimes (ii (length statelist))
      (setf (gethash (nth ii statelist) ht) ii))
    ht))
(defun symbolize (atree)
;  (format t "symbolize: atree=~A~%" atree)
  (cond ((definition-p atree) 
         (list (assym1 (show (definition-name atree)))
               (symbolize (definition-sform atree))))
        ((transition-p atree)
         (list (assym1 (show (transition-name atree)))
               (symbolize (transition-value atree))))
        ((gname-p atree) (assym1 (show atree)))
        ((atom atree) atree)
        ((consp atree) (mapcar #'symbolize atree))
        (t (error "symbolize: error")))) 
@

<<>>=
(defun compile-transitions (translist statelist)
  (let* ((dlist (dup NIL (length statelist)))
         (frst NIL)
         (tost NIL)
         (ltable (make-hash-table :test #'equal)))
    (dotimes (ii (length statelist))
      (setf (gethash (nth ii statelist) ltable) ii))
    (dolist (tr translist) 
      (setf frst (gethash (transition-from tr) ltable))
      (setf tost (gethash (transition-to tr) ltable))
      (if (and (null frst) (not (source-p (transition-from tr))))
          (progn
           (format t "transition=~A~%" tr)
           (format t "frst=~A; (source-p frst)=~A~%" frst (source-p frst))
           (error "from state not in statelist")))
      (if (and (null tost) (not (sink-p (transition-to tr))))
          (progn
           (format t "transition=~A~%" tr)
           (error "to-state not in statelist")))
      (if (not (source-p (transition-from tr)))
          (setf (nth frst dlist) 
                (cons (cons '- (list (symbolize (transition-name tr))))
                      (nth frst dlist)))) 
      (if (not (sink-p (transition-to tr)))
          (setf (nth tost dlist)
                (cons (transition-name tr) (nth tost dlist)))))
    (dotimes (ii (length dlist))
      (setf (nth ii dlist) (cons '+ (nth ii dlist))))
    (cons 'vector dlist)))
@

<<>>=
(defun has-factor (expression state)
  (find state expression))
(defun factor-out (expression state)
  (cond ((eq (car expression) '*)
         (let* ((e1 (remove state expression :count 1)))
           (cond ((eq (length e1) 2)
                  (cadr e1))
                 (t e1))))
        (t (error "factor-out: not a product"))))
(defun compile-transitions-factored (translist statelist)
  (let* ((positive (dup NIL (length statelist)))
         (negative (dup NIL (length statelist)))
         (dlist (dup NIL (length statelist)))
         (negfactors (dup NIL (length statelist)))
         (frst NIL)
         (tost NIL)
         (cursym NIL)
         (nf NIL)
         (sta NIL)
         (fc NIL)
         (rmu NIL)
         (expr NIL)
         (ans NIL)
         (rans NIL)
         (nans NIL)
         (maxargs NIL)
         (tmpname NIL)
         (ltable (make-hash-table :test #'equal)))
    (dotimes (ii (length statelist))
      (setf (gethash (nth ii statelist) ltable) ii))
    (dolist (tr translist) 
      (setf frst (gethash (transition-from tr) ltable))
      (setf tost (gethash (transition-to tr) ltable))
      (if (not (source-p (transition-from tr)))
          (cond ((has-factor (transition-value tr) (transition-from tr))
                 (progn
                  (setf cursym (transition-name tr))
                  (setf (get cursym 'factor)  
                        (factor-out (transition-value tr) 
                                    (transition-from tr)))
                  (setf (nth frst negfactors) 
                        (cons cursym (nth frst negfactors)))
                  (setf (nth frst dlist) 
                        (cons (cons '- (list (transition-name tr)))
                              (nth frst dlist)))))
                (t (setf (nth frst negative)
                         (cons (transition-name tr) 
                               (nth frst negative))))))
      (if (not (sink-p (transition-to tr)))
          (setf (nth tost dlist)
                (cons (transition-name tr)
                      (nth tost dlist)))))
    (dotimes (ii (length negfactors))
      (setf nf (nth ii negfactors))
      (setf sta (nth ii statelist))
      (setf fc (mapcar #'(lambda (x) (get x 'factor)) nf))
      (setf rmu (if (null fc) NIL (list 'RTRANS ':size sta ':prob fc ':deltat 'DELTAT :npars (length fc))))
      (setf expr (list nf rmu))
      (setf ans (cons expr ans)))
    (setf ans (remove-if #'(lambda (x) (null (car x))) ans))
    (setf rans (reverse (copy-list ans)))
    (setf *rans* rans)
    (setf nans NIL)
    (setf maxargs 0)
    (dolist (rr rans)
      (progn
       (setf tmpname (new-gname (gensym)))
       (setf (get tmpname 'pointer) T)
       (if (> (length (car rr)) maxargs)
           (setf maxargs (length (car rr))))
       (dotimes (ii (length (car rr)))
         (setf nans (cons (list (nth ii (car rr))
                                (list 'AREF tmpname ii))
                          nans)))
       (setf nans (cons (list tmpname (cadr rr)) nans))))
    (setf nans (cons (list 'speci@l-mark maxargs) nans))
    (dotimes (ii (length dlist))
      (setf (nth ii dlist) (cons '+ (nth ii dlist))))
    (list nans (cons 'vector dlist))))
@

<<>>=
(defun to-lisp-ode (amodel)
  (let ((pars (pa amodel))
        (states (st amodel))
        (dtlist (remove-if-not #'dtfilter (topo-sort amodel)))
        (zzzk (prin1 "done with topo-sort")))
    (let ((def1 (append (to-def 'PARAM pars) 
                        (to-def 'STATE states)
                        (symbolize dtlist))))
      (list 'LAMBDA 
            (list 'PARAM 'STATE 'SIMULATION-TIME) 
            (list 'LET* def1 (compile-transitions 
                              (model-transitions amodel)
                              (st amodel)))))))
(defun to-lisp-sde (amodel)
  (let ((pars (pa amodel))
        (states (st amodel))
        (dtlist (remove-if-not #'dtfilter (topo-sort amodel))))
    (let* ((def1 (append (to-def 'PARAM pars) 
                         (to-def 'STATE states)
                         (symbolize dtlist)))
           (ctf (compile-transitions-factored (model-transitions amodel) (st amodel)))
           (def2 (car ctf))
           (def3 (append def1 def2)))
      (list 'LAMBDA 
            (list 'PARAM 'STATE 'SIMULATION-TIME 'DELTAT) 
            (list 'LET* def3 (cadr ctf))))))
@
One way to use \verb+to-lisp-ode+:
\begin{verbatim}
(defvar z (to-lisp-ode mod1))
(eval (list 'function z))
\end{verbatim}

\section*{Compilation to C}
The functions that follow compile lambda expressions featuring a minute subset
of Lisp into C.
<<>>=
(defun ^ (aa bb) (exp (* bb (log aa))))
@

<<>>=
(defun iseq (aa bb)
  (do ((ii bb (- ii 1))
       (ans () (cons ii ans)))
      ((< ii aa) ans)))       
@

<<>>=
(defun s^ (aa bb)
  (cond ((equal bb 0) 1)
        ((equal bb 1) aa)
        ((and (numberp aa) (numberp bb)) (^ aa bb))
        (t (list '^ aa bb))))
@

<<>>=
;; split long lists before calling s*
(defun ss* (alist)
  (cond ((atom alist) alist)
        ((equal (length alist) 1) (car alist))
        (t (let* ((thelists (splitdown alist 10))
                  (ss (mapcar 
                       #'(lambda (onelist) 
                          (apply #'s* onelist)) 
                       thelists)))
             (funcall #'ss* ss)))))
@

<<>>=
(defun s* (&rest arglist)
  (let ((len (length arglist)))
    (cond ((> len 11)
           (let ((splits (split arglist (round (/ len 2)))))
             (s* (apply 's* (car splits))
                 (apply 's* (cadr splits)))))
          (t (gentimes-aux arglist () 1)))))
@

<<>>=
(defun sdot (vec1 vec2)
  (apply #'s+ (mapcar #'s* vec1 vec2)))
@

<<>>=
(defun gentimes-aux (ins syms numbers)
  (cond ((null ins)
         (cond ((null syms) numbers)
               ((equal numbers 1) 
                (cond ((null (cdr syms)) (car syms))
                      (t (cons '* syms))))
               (t (cons '* (cons numbers syms)))))
        ((equal (car ins) 0) 0)
        ((equal (car ins) 1) 
         (gentimes-aux (cdr ins) syms numbers))
        ((numberp (car ins))
         (gentimes-aux (cdr ins) syms (* numbers (car ins))))
        (t
         (gentimes-aux (cdr ins) (cons (car ins) syms) numbers))))
@

This splits list into two parts, first nn elements to first
output list.  
<<>>=
(defun split (alist nn)
  (do ((curlist alist (cdr curlist))
       (ii 0 (+ 1 ii))
       (ans () (cons (car curlist) ans)))
      ((>= ii nn) (list (reverse ans) curlist))))
@

Split alist into a list of lists each no larger than nn elements
<<>>=
(defun splitdown (alist nn)
  (do* ((tsplit () (split curlist nn))
        (curlist alist (cadr tsplit))
        (len (length alist) (- len nn))
        (accum () (cons (car tsplit) accum)))
       ((<= len nn) (reverse (cons curlist accum)))))
@

Split long lists before calling s+:
<<>>=
(defun ss+ (alist)
  (cond ((atom alist) alist)
        ((equal (length alist) 1) (car alist))
        (t (let* ((thelists (splitdown alist 10))
                  (ss (mapcar 
                       #'(lambda (onelist) 
                          (apply #'s+ onelist)) 
                       thelists)))
             (funcall #'ss+ ss)))))
@

<<>>=
(defun s+ (&rest arglist)
  (let ((len (length arglist)))
    (cond ((> len 31)
           (let ((splits (split arglist (round (/ len 2)))))
             (s+ (apply 's+ (car splits))
                 (apply 's+ (cadr splits)))))
          (t (genplus-aux arglist () 0)))))
@

<<>>=
(defun genplus-aux (ins syms numbers)
  (cond ((null ins)
         (cond ((null syms) numbers)
               ((equal numbers 0) 
                (cond ((null (cdr syms)) (car syms))
                      (t (cons '+ syms))))
               (t (cons '+ (cons numbers syms)))))
        ((equal (car ins) 0) 
         (genplus-aux (cdr ins) syms numbers))
        ((numberp (car ins))
         (genplus-aux (cdr ins) syms (+ numbers (car ins))))
        (t
         (genplus-aux (cdr ins) (cons (car ins) syms) numbers))))
@

<<>>=
(defun s- (aa bb)
  (cond ((and (numberp aa) (numberp bb)) (- aa bb))
        ((eql bb 0) aa)
        ((eql aa 0) (list '- bb))
        (t (list '- aa bb))))
@

<<>>=
(defun s-unary-minus (aa)
  (cond ((numberp aa) (- aa))
        ((equal aa 0) 0)
        (t (list '- aa))))
@

<<>>=
(defun s/ (aa bb)
  (cond ((equal bb 0) (error "zerodivide"))
        ((and (numberp aa) (numberp bb)) (/ aa bb))
        (t (list '/ aa bb))))
@

<<>>=
(defun s-eval-log (expr)
  (cond ((numberp expr) (log expr))
        ((symbolp expr) (list 'log expr))
        ((equal (car expr) '*)
         (cons '+ (mapcar #'s-eval-log (cdr expr))))
        ((equal (car expr) '^)
         (list '* (caddr expr) (s-eval-log (cadr expr))))
        (t (list 'log expr))))
@

<<>>=
(defstruct cfunction
  hdr
  dclrs
  stmts
  closer
)
(defstruct cvar
  name
  type
  dimension
)
@
<<>>=
(defun sh (x) (progn (format t "sh: ~A~%" x) x))
@
<<>>=
(defun all-keys (ht)
  (let ((ans ()))
    (maphash #'(lambda (x y) (setf ans (cons x ans))) ht)
    ans))
@

<<>>=
(defun nodashorat (ast)
  (let* ((str (remove "["
                (remove ","
                 (remove "]"
                  (remove "_"
                   (remove "-"
                    (remove "@" (string (symbolize ast)) :test #'string-equal)
                    :test #'string-equal)
                   :test #'string-equal)
                  :test #'string-equal)
                 :test #'string-equal)
                :test #'string-equal)))
    str))
@

<<>>=
(defun compile-to-c (amodel &key lisp-form)
  (let* ((vs (car (cdaddr lisp-form)))
         (tmp0 (setf *nr-tmps-c* 0))
         (osyms NIL)
         (lisp-form-copy lisp-form)
         (bodyvars (mapcar #'car vs))
         (markloc (position-if #'(lambda (x) (eq x 'speci@l-mark)) bodyvars))
         (othervars (if (null markloc) nil (cdr (nthcdr markloc bodyvars))))
         (vecsize (if (null markloc) nil (cadr (nth markloc (car (cdaddr lisp-form))))))
         (statelist (st amodel))
         (parlist (pa amodel))
         (trlist (mapcar #'(lambda (x) (symbolize (transition-name x))) (model-transitions amodel)))
         (deflist (mapcar #'(lambda (x) (symbolize (definition-name x))) (model-definitions amodel))))
    (dribble "csysode.1")
    (format t "#include <stdio.h>~%")
    (format t "#include <math.h>~%")
    (format t "#include csystmp.h~%") ; !!!Fix
    (format t "#define LOG\(x\) log\(\(double\) x\)~%")
    ; now, add preprocessor substitutions for params and states
    (do ((rec (cdr vs) (cdr rec))
         (top (format t "#define ~A ~A[~A]~%" (nodashorat (caar vs))
                         (cadr (cadar vs)) 
                         (caddr (cadar vs))) 
              (format t "#define ~A ~A[~A]~%" (nodashorat (caar rec))
                         (cadr (cadar rec)) 
                         (caddr (cadar rec)))))
      ((not (or (eq (cadr (cadar rec)) 'param) (eq (cadr (cadar rec)) 'state))) (setf osyms rec)))
    (dribble "csysode.1")
    (format t "void csysode(double *STATE, double *PARAM, double *SIMTIME, double *DELTAT, double *dxdt){~%")
    (format t "/*~%")
    (format t "  double STATE[~A] ;~%" (length (st amodel)))
    (format t "  double PARAM[~A] ;~%" (length (pa amodel)))
    (format t "  double dxdt[~A] ;~%" (length (st amodel)))
    (format t "*/~%")
    (format t "  double tmp[CSYSTMP] ;~%")
    ; make dictionary of names to numbers; state or param names not included
    (setf symtab (make-hash-table :test #'equal))
    (setf curno 0)
    (dolist (tr trlist)
      ; (format t "~A ~A~%" (nodashorat tr) curno)
      (setf (gethash (nodashorat tr) symtab) curno)
      (incf curno)) 
    (dolist (df deflist)
      ; (format t "~A ~A~%" (nodashorat df) curno)
      (setf (gethash (nodashorat df) symtab) curno)
      (incf curno))
    (setf vecsym "u")
    (format t "  double ~A[~A];~%" vecsym curno)
    ; modify: leave out state/param clauses in lisp-form
    (setf (car (cdaddr lisp-form-copy)) osyms)
    (infix lisp-form-copy vecsym symtab)
    ;(infix lisp-form-copy NIL symtab)
    (format t "}~%")
    (dribble)
    (dribble "csystmp.h")
    (format t "#define CSYSTMP ~A~%" *nr-tmps-c*) 
    (dribble)
    (dribble "sysaux.r")
    ; COMMENT: do NOT change the next line
    (format t "/* R sysode.r function ******~%")
    (format t "sysode <- function(t,y,params) {~%")
    (format t "  list(.C(\"csysode\", as.double(y), as.double(params),~%")
    (format t "          as.double(t), as.double(0.0),~%")
    (format t "          ans=as.double(rep(0,~A)))$ans,~%" (length (st amodel)))
    (format t "       NULL)~%}~%")

    ; COMMENT: do NOT change the next line
    (format t "************** end sysode.r ********/~%")
    ; COMMENT: do NOT change the next line
    (format t "/* R dosim.r function ******~%")
    (format t "dyn\.load\(\"./csysode.so\"\);~%")
    (format t "require\(\"odesolve\"\);~%")
    (format t "source\(\"./sysode.r\"\);~%") 
    (format t "# Parameter value section~%") 
    (dotimes (ii (length parlist))
      (format t "~A <- ???;~%" (nodashorat (nth ii parlist))))
    (format t "~%# Initial value section~%")
    (dotimes (ii (length statelist))
      (format t "~A~A <- ???;~%" (nodashorat (nth ii statelist)) "in"))
    (format t "# End time~%")
    (format t "tend <- ???;~%")
    (format t "~%")
    (format t "param <- rep\(0,~A\);~%" (length parlist))
    (dotimes (ii (length parlist))
      (format t "param\[~A\] <- ~A;~%" (+ 1 ii) (nodashorat (nth ii parlist)))) 
    (format t "yin <- rep\(0,~A\);~%" (length statelist))
    (dotimes (ii (length statelist))
      (format t "yin\[~A\] <- ~Ain;~%" (+ 1 ii) (nodashorat (nth ii statelist))))
    (format t "tt<-seq\(0,tend,by=0.125\)~%")
    (format t "yout <- lsoda\(yin,tt,sysode,param,rtol=1e-4,atol=1e-3\);~%")
    (format t "# states:~%") 
    (dotimes (ii (length statelist))
      (format t "~A <- yout\[,~A\];~%" (nodashorat (nth ii statelist)) (+ 2 ii)))
    ; COMMENT: do NOT change the next line either
    (format t "************** end dosim.r ********/~%")

    (format t "/* R boilerplate; can be cpp~%")
    (dotimes (ii (length parlist))
      (format t "params\[~A\] <- ~A;~%" (+ 1 ii) (nodashorat (nth ii parlist)))) 
    (format t "*/~%")
    (format t "~%/* R boilerplate; can be cpp~%")
    (dotimes (ii (length parlist))
      (format t "~A <- params\[~A\];~%" (nodashorat (nth ii parlist)) (+ 1 ii)))
    (format t "*/~%")
    (format t "~%/* R boilerplate; can be cpp~%")
    (dotimes (ii (length statelist))
      (format t "yy\[~A\] <- ~A;~%" (+ 1 ii) (nodashorat (nth ii statelist))))
    (format t "*/~%")
    (format t "~%/* R boilerplate; can be cpp~%")
    (dotimes (ii (length statelist))
      (format t "~A <- yy\[~A\];~%" (nodashorat (nth ii statelist)) (+ 1 ii)))
    (format t "*/~%")
    (dribble)
))
@

<<>>=
(defun infix (expr &optional (sym NIL) (htab NIL))
  (cond ((numberp expr)
	 (format t "~A" expr))
        ((not (listp expr)) 
	 (format t "(")
         (if (null sym)
             (format t "~A" (nodashorat expr))
           (let ((onum (gethash (nodashorat expr) htab)))
             (if (null onum)
                 (format t "~A" (nodashorat expr))
               (format t "~A[~A]" sym onum))))
	 (format t ")"))
	((equalp (car expr) 'aref)
	 (progn
	   (format t "(")
	   (format t "~A" (nodashorat (second expr)))
	   (format t "[")
	   (format t "~A" (third expr))
	   (format t "])")))
        ((equalp (car expr) '+)
         (if (<= (length (cdr expr)) 10)  
  	     (progn
	      (format t "(")
	      (infix (cadr expr) sym htab)
	      (dolist (cur (cddr expr))
	        (format t "+")
	        (infix cur sym htab))
	      (format t ")"))
           (progn
            (format t "(")
            (format t "tmp[~A]=0.0," *nr-tmps-c*)        
            (dolist (cur (cdr expr))
              (format t "tmp[~A]+=" *nr-tmps-c*)
              (infix cur sym htab)
              (format t ","))
            (format t "tmp[~A])~%" *nr-tmps-c*)
            (incf *nr-tmps-c*))))
        ((equalp (car expr) '*)
	 (progn
	   (format t "(")
	   (infix (cadr expr) sym htab)
	   (dolist (cur (cddr expr))
	     (format t "*")
	     (infix cur sym htab))
	   (format t ")")))
        ((equalp (car expr) '-)
	 (cond ((equalp (length expr) 2)  ; unary minus
 		(format t "(-")
		(infix (cadr expr) sym htab)
		(format t ")"))
	       (t (progn
	           (format t "(")
	           (infix (cadr expr) sym htab)
	           (dolist (cur (cddr expr))
	            (format t "-")
	            (infix cur sym htab)
	            (format t ")"))))))
;        ((equalp (car expr) '/)
;	 (progn
;	   (format t "(")
;	   (infix (cadr expr) sym htab)
;	   (dolist (cur (cddr expr))
;	     (format t "/")
;	     (infix cur sym htab))
;	   (format t ")")))
        ((equalp (car expr) '/)
	 (progn
	   (format t "((")
	   (infix (cadr expr) sym htab)
	   (format t ")")
	   (format t "/")
	   (format t "(")
	   (infix (caddr expr) sym htab)
	   (format t "))")))
	((equalp (car expr) 'the)
	 (progn
	   (infix (third expr) sym htab)))
	((equalp (car expr) 'lambda)
	 (progn
	   (dolist (cur (cddr expr))
	     (cond ((equalp (car cur) 'declare) (infix cur sym htab))
		   (t (infix cur sym htab))))))
	((equalp (car expr) 'vector)
	 (let ((ii 0))
	   (format t "dxdt[~A]=" ii)
           (incf ii)
	   (infix (cadr expr) sym htab)
	   (dolist (cur (cddr expr))
	     (format t " ;") (terpri)
             (format t "dxdt[~A]=" ii)
             (incf ii)
	     (infix cur sym htab))
	   (format t " ;") (terpri)))
	((equalp (car expr) 'let*)	; caution--this won't work if nested
	 (progn
	   (dolist (defn (cadr expr))
             (if (not (eq (car defn) 'speci@l-mark))
                 (progn
                   (if (null sym)
  	               (format t "~A=" (nodashorat (car defn)))
                     (let ((onum (gethash (nodashorat (car defn)) htab)))
                       (if (null onum)
                           (format t "~A=" (nodashorat (car defn)))
                         (format t "~A[~A]=" sym onum))))
	           (infix (cadr defn) sym htab)
	           (format t ";~%"))))
	   (infix (third expr) sym htab)))
        ((equalp (car expr) 'RTRANS)
         (progn
           (format t "(")
           (format t "rtrans(outarg,(")
           (dotimes (ii (length (fifth expr)))
             (format t "inarg[~A]=~A," ii (nth ii (fifth expr))))
           (format t "inarg),~A,~A,~A))" (third expr) (seventh expr) (ninth expr))))
	(t
         (progn
   	   (format t "(")
           (if (null sym)
               (format t "~A" (nodashorat (car expr)))
             (let ((onum (gethash (nodashorat (car expr)) htab)))
               (if (null onum) 
                   (format t "~A" (nodashorat (car expr)))
                 (format t "~A[~A]" sym onum))))
	   (format t "(")
	   (infix (cadr expr) sym htab)
	   (dolist (cur (cddr expr))
	     (format t ",")
	     (infix cur sym htab))
	   (format t "))")))))
@

Finally:
<<>>=
(defun sha (ht) (maphash #'(lambda (x y) (print x)(print y)(print "--")) ht))

(defmacro finish (amodel)
  (let ((zu (gensym)))
    `(progn
      (defvar ,zu (to-lisp-ode ,amodel))
      (compile-to-c ,amodel :lisp-form ,zu))))
@

\newpage
\section*{Tuberculosis model}

Here, we will build the TB model.  
\begin{verbatim}
;    __ABCD-E means the number of people whose wild type
;     indicator status is A, whose INH resistant 
;       subpopulation indicator status is B, whose RIF subpopulation indicator
;       status is C, and whose MDR subpopulation indicator status is D, and
;       whose infection is dominated by E.  A %in% c(0,1), B %in% c(0,1),
;       C %in% c(0,1), D %in% c(0,1); E %in% 1:4
\end{verbatim}
<<>>=
(defvar tbmodel (new-model 'tbmodel))
@

<<>>=
(defvar new-fast-latent-states 
  (list
    (new-gname 'LL-1000-1)
    (new-gname 'LL-0100-2)
    (new-gname 'LL-0010-3)
    (new-gname 'LL-0001-4)
    (new-gname 'LL-1100-1)
    (new-gname 'LL-1100-2)
    (new-gname 'LL-1010-1)
    (new-gname 'LL-1010-3)
    (new-gname 'LL-1001-1)
    (new-gname 'LL-1001-4)
    (new-gname 'LL-0110-2)
    (new-gname 'LL-0110-3)
    (new-gname 'LL-0101-2)
    (new-gname 'LL-0101-4)
    (new-gname 'LL-0011-3)
    (new-gname 'LL-0011-4)
    (new-gname 'LL-0111-2)
    (new-gname 'LL-0111-3)
    (new-gname 'LL-0111-4)
    (new-gname 'LL-1011-1)
    (new-gname 'LL-1011-3)
    (new-gname 'LL-1011-4)
    (new-gname 'LL-1101-1)
    (new-gname 'LL-1101-2)
    (new-gname 'LL-1101-4)
    (new-gname 'LL-1110-1)
    (new-gname 'LL-1110-2)
    (new-gname 'LL-1110-3)
    (new-gname 'LL-1111-1)
    (new-gname 'LL-1111-2)
    (new-gname 'LL-1111-3)
    (new-gname 'LL-1111-4)))
@

<<>>=
(defvar new-slow-latent-states 
  (list
    (new-gname 'MM-1000-1)
    (new-gname 'MM-0100-2)
    (new-gname 'MM-0010-3)
    (new-gname 'MM-0001-4)
    (new-gname 'MM-1100-1)
    (new-gname 'MM-1100-2)
    (new-gname 'MM-1010-1)
    (new-gname 'MM-1010-3)
    (new-gname 'MM-1001-1)
    (new-gname 'MM-1001-4)
    (new-gname 'MM-0110-2)
    (new-gname 'MM-0110-3)
    (new-gname 'MM-0101-2)
    (new-gname 'MM-0101-4)
    (new-gname 'MM-0011-3)
    (new-gname 'MM-0011-4)
    (new-gname 'MM-0111-2)
    (new-gname 'MM-0111-3)
    (new-gname 'MM-0111-4)
    (new-gname 'MM-1011-1)
    (new-gname 'MM-1011-3)
    (new-gname 'MM-1011-4)
    (new-gname 'MM-1101-1)
    (new-gname 'MM-1101-2)
    (new-gname 'MM-1101-4)
    (new-gname 'MM-1110-1)
    (new-gname 'MM-1110-2)
    (new-gname 'MM-1110-3)
    (new-gname 'MM-1111-1)
    (new-gname 'MM-1111-2)
    (new-gname 'MM-1111-3)
    (new-gname 'MM-1111-4)))
@

<<>>=
(defvar smearnegative-states
  (list
    (new-gname 'SS-1000-1)
    (new-gname 'SS-0100-2)
    (new-gname 'SS-0010-3)
    (new-gname 'SS-0001-4)
    (new-gname 'SS-1100-1)
    (new-gname 'SS-1100-2)
    (new-gname 'SS-1010-1)
    (new-gname 'SS-1010-3)
    (new-gname 'SS-1001-1)
    (new-gname 'SS-1001-4)
    (new-gname 'SS-0110-2)
    (new-gname 'SS-0110-3)
    (new-gname 'SS-0101-2)
    (new-gname 'SS-0101-4)
    (new-gname 'SS-0011-3)
    (new-gname 'SS-0011-4)
    (new-gname 'SS-0111-2)
    (new-gname 'SS-0111-3)
    (new-gname 'SS-0111-4)
    (new-gname 'SS-1011-1)
    (new-gname 'SS-1011-3)
    (new-gname 'SS-1011-4)
    (new-gname 'SS-1101-1)
    (new-gname 'SS-1101-2)
    (new-gname 'SS-1101-4)
    (new-gname 'SS-1110-1)
    (new-gname 'SS-1110-2)
    (new-gname 'SS-1110-3)
    (new-gname 'SS-1111-1)
    (new-gname 'SS-1111-2)
    (new-gname 'SS-1111-3)
    (new-gname 'SS-1111-4)))
@

<<>>=
(defvar smearpositive-states
  (list
    (new-gname 'PP-1000-1)
    (new-gname 'PP-0100-2)
    (new-gname 'PP-0010-3)
    (new-gname 'PP-0001-4)
    (new-gname 'PP-1100-1)
    (new-gname 'PP-1100-2)
    (new-gname 'PP-1010-1)
    (new-gname 'PP-1010-3)
    (new-gname 'PP-1001-1)
    (new-gname 'PP-1001-4)
    (new-gname 'PP-0110-2)
    (new-gname 'PP-0110-3)
    (new-gname 'PP-0101-2)
    (new-gname 'PP-0101-4)
    (new-gname 'PP-0011-3)
    (new-gname 'PP-0011-4)
    (new-gname 'PP-0111-2)
    (new-gname 'PP-0111-3)
    (new-gname 'PP-0111-4)
    (new-gname 'PP-1011-1)
    (new-gname 'PP-1011-3)
    (new-gname 'PP-1011-4)
    (new-gname 'PP-1101-1)
    (new-gname 'PP-1101-2)
    (new-gname 'PP-1101-4)
    (new-gname 'PP-1110-1)
    (new-gname 'PP-1110-2)
    (new-gname 'PP-1110-3)
    (new-gname 'PP-1111-1)
    (new-gname 'PP-1111-2)
    (new-gname 'PP-1111-3)
    (new-gname 'PP-1111-4)))
  
(defvar old-fast-latent-states
  (list
    (new-gname 'LLPRIME-1000-1)
    (new-gname 'LLPRIME-0100-2)
    (new-gname 'LLPRIME-0010-3)
    (new-gname 'LLPRIME-0001-4)
    (new-gname 'LLPRIME-1100-1)
    (new-gname 'LLPRIME-1100-2)
    (new-gname 'LLPRIME-1010-1)
    (new-gname 'LLPRIME-1010-3)
    (new-gname 'LLPRIME-1001-1)
    (new-gname 'LLPRIME-1001-4)
    (new-gname 'LLPRIME-0110-2)
    (new-gname 'LLPRIME-0110-3)
    (new-gname 'LLPRIME-0101-2)
    (new-gname 'LLPRIME-0101-4)
    (new-gname 'LLPRIME-0011-3)
    (new-gname 'LLPRIME-0011-4)
    (new-gname 'LLPRIME-0111-2)
    (new-gname 'LLPRIME-0111-3)
    (new-gname 'LLPRIME-0111-4)
    (new-gname 'LLPRIME-1011-1)
    (new-gname 'LLPRIME-1011-3)
    (new-gname 'LLPRIME-1011-4)
    (new-gname 'LLPRIME-1101-1)
    (new-gname 'LLPRIME-1101-2)
    (new-gname 'LLPRIME-1101-4)
    (new-gname 'LLPRIME-1110-1)
    (new-gname 'LLPRIME-1110-2)
    (new-gname 'LLPRIME-1110-3)
    (new-gname 'LLPRIME-1111-1)
    (new-gname 'LLPRIME-1111-2)
    (new-gname 'LLPRIME-1111-3)
    (new-gname 'LLPRIME-1111-4)))
@

<<>>=
(defvar old-slow-latent-states
  (list
    (new-gname 'MMPRIME-1000-1)
    (new-gname 'MMPRIME-0100-2)
    (new-gname 'MMPRIME-0010-3)
    (new-gname 'MMPRIME-0001-4)
    (new-gname 'MMPRIME-1100-1)
    (new-gname 'MMPRIME-1100-2)
    (new-gname 'MMPRIME-1010-1)
    (new-gname 'MMPRIME-1010-3)
    (new-gname 'MMPRIME-1001-1)
    (new-gname 'MMPRIME-1001-4)
    (new-gname 'MMPRIME-0110-2)
    (new-gname 'MMPRIME-0110-3)
    (new-gname 'MMPRIME-0101-2)
    (new-gname 'MMPRIME-0101-4)
    (new-gname 'MMPRIME-0011-3)
    (new-gname 'MMPRIME-0011-4)
    (new-gname 'MMPRIME-0111-2)
    (new-gname 'MMPRIME-0111-3)
    (new-gname 'MMPRIME-0111-4)
    (new-gname 'MMPRIME-1011-1)
    (new-gname 'MMPRIME-1011-3)
    (new-gname 'MMPRIME-1011-4)
    (new-gname 'MMPRIME-1101-1)
    (new-gname 'MMPRIME-1101-2)
    (new-gname 'MMPRIME-1101-4)
    (new-gname 'MMPRIME-1110-1)
    (new-gname 'MMPRIME-1110-2)
    (new-gname 'MMPRIME-1110-3)
    (new-gname 'MMPRIME-1111-1)
    (new-gname 'MMPRIME-1111-2)
    (new-gname 'MMPRIME-1111-3)
    (new-gname 'MMPRIME-1111-4)))
  
(defvar incomplete-cure-states
  (list
    (new-gname 'RR-1000-1)
    (new-gname 'RR-0100-2)
    (new-gname 'RR-0010-3)
    (new-gname 'RR-0001-4)
    (new-gname 'RR-1100-1)
    (new-gname 'RR-1100-2)
    (new-gname 'RR-1010-1)
    (new-gname 'RR-1010-3)
    (new-gname 'RR-1001-1)
    (new-gname 'RR-1001-4)
    (new-gname 'RR-0110-2)
    (new-gname 'RR-0110-3)
    (new-gname 'RR-0101-2)
    (new-gname 'RR-0101-4)
    (new-gname 'RR-0011-3)
    (new-gname 'RR-0011-4)
    (new-gname 'RR-0111-2)
    (new-gname 'RR-0111-3)
    (new-gname 'RR-0111-4)
    (new-gname 'RR-1011-1)
    (new-gname 'RR-1011-3)
    (new-gname 'RR-1011-4)
    (new-gname 'RR-1101-1)
    (new-gname 'RR-1101-2)
    (new-gname 'RR-1101-4)
    (new-gname 'RR-1110-1)
    (new-gname 'RR-1110-2)
    (new-gname 'RR-1110-3)
    (new-gname 'RR-1111-1)
    (new-gname 'RR-1111-2)
    (new-gname 'RR-1111-3)
    (new-gname 'RR-1111-4)))
@

<<>>=
(defvar treated-smneg-states-standard
  (list
    (new-gname 'YYNEG-1000-1)
    (new-gname 'YYNEG-0100-2)
    (new-gname 'YYNEG-0010-3)
    (new-gname 'YYNEG-0001-4)
    (new-gname 'YYNEG-1100-1)
    (new-gname 'YYNEG-1100-2)
    (new-gname 'YYNEG-1010-1)
    (new-gname 'YYNEG-1010-3)
    (new-gname 'YYNEG-1001-1)
    (new-gname 'YYNEG-1001-4)
    (new-gname 'YYNEG-0110-2)
    (new-gname 'YYNEG-0110-3)
    (new-gname 'YYNEG-0101-2)
    (new-gname 'YYNEG-0101-4)
    (new-gname 'YYNEG-0011-3)
    (new-gname 'YYNEG-0011-4)
    (new-gname 'YYNEG-0111-2)
    (new-gname 'YYNEG-0111-3)
    (new-gname 'YYNEG-0111-4)
    (new-gname 'YYNEG-1011-1)
    (new-gname 'YYNEG-1011-3)
    (new-gname 'YYNEG-1011-4)
    (new-gname 'YYNEG-1101-1)
    (new-gname 'YYNEG-1101-2)
    (new-gname 'YYNEG-1101-4)
    (new-gname 'YYNEG-1110-1)
    (new-gname 'YYNEG-1110-2)
    (new-gname 'YYNEG-1110-3)
    (new-gname 'YYNEG-1111-1)
    (new-gname 'YYNEG-1111-2)
    (new-gname 'YYNEG-1111-3)
    (new-gname 'YYNEG-1111-4)))
@

<<>>=
(defvar treated-smneg-states-specialized
  (list
    (new-gname 'WWNEG-1000-1)
    (new-gname 'WWNEG-0100-2)
    (new-gname 'WWNEG-0010-3)
    (new-gname 'WWNEG-0001-4)
    (new-gname 'WWNEG-1100-1)
    (new-gname 'WWNEG-1100-2)
    (new-gname 'WWNEG-1010-1)
    (new-gname 'WWNEG-1010-3)
    (new-gname 'WWNEG-1001-1)
    (new-gname 'WWNEG-1001-4)
    (new-gname 'WWNEG-0110-2)
    (new-gname 'WWNEG-0110-3)
    (new-gname 'WWNEG-0101-2)
    (new-gname 'WWNEG-0101-4)
    (new-gname 'WWNEG-0011-3)
    (new-gname 'WWNEG-0011-4)
    (new-gname 'WWNEG-0111-2)
    (new-gname 'WWNEG-0111-3)
    (new-gname 'WWNEG-0111-4)
    (new-gname 'WWNEG-1011-1)
    (new-gname 'WWNEG-1011-3)
    (new-gname 'WWNEG-1011-4)
    (new-gname 'WWNEG-1101-1)
    (new-gname 'WWNEG-1101-2)
    (new-gname 'WWNEG-1101-4)
    (new-gname 'WWNEG-1110-1)
    (new-gname 'WWNEG-1110-2)
    (new-gname 'WWNEG-1110-3)
    (new-gname 'WWNEG-1111-1)
    (new-gname 'WWNEG-1111-2)
    (new-gname 'WWNEG-1111-3)
    (new-gname 'WWNEG-1111-4)))
@

<<>>=
(defvar treated-smpos-states-standard
  (list
    (new-gname 'YYPOS-1000-1)
    (new-gname 'YYPOS-0100-2)
    (new-gname 'YYPOS-0010-3)
    (new-gname 'YYPOS-0001-4)
    (new-gname 'YYPOS-1100-1)
    (new-gname 'YYPOS-1100-2)
    (new-gname 'YYPOS-1010-1)
    (new-gname 'YYPOS-1010-3)
    (new-gname 'YYPOS-1001-1)
    (new-gname 'YYPOS-1001-4)
    (new-gname 'YYPOS-0110-2)
    (new-gname 'YYPOS-0110-3)
    (new-gname 'YYPOS-0101-2)
    (new-gname 'YYPOS-0101-4)
    (new-gname 'YYPOS-0011-3)
    (new-gname 'YYPOS-0011-4)
    (new-gname 'YYPOS-0111-2)
    (new-gname 'YYPOS-0111-3)
    (new-gname 'YYPOS-0111-4)
    (new-gname 'YYPOS-1011-1)
    (new-gname 'YYPOS-1011-3)
    (new-gname 'YYPOS-1011-4)
    (new-gname 'YYPOS-1101-1)
    (new-gname 'YYPOS-1101-2)
    (new-gname 'YYPOS-1101-4)
    (new-gname 'YYPOS-1110-1)
    (new-gname 'YYPOS-1110-2)
    (new-gname 'YYPOS-1110-3)
    (new-gname 'YYPOS-1111-1)
    (new-gname 'YYPOS-1111-2)
    (new-gname 'YYPOS-1111-3)
    (new-gname 'YYPOS-1111-4)))
@

<<>>=
(defvar treated-smpos-states-specialized
  (list
    (new-gname 'WWPOS-1000-1)
    (new-gname 'WWPOS-0100-2)
    (new-gname 'WWPOS-0010-3)
    (new-gname 'WWPOS-0001-4)
    (new-gname 'WWPOS-1100-1)
    (new-gname 'WWPOS-1100-2)
    (new-gname 'WWPOS-1010-1)
    (new-gname 'WWPOS-1010-3)
    (new-gname 'WWPOS-1001-1)
    (new-gname 'WWPOS-1001-4)
    (new-gname 'WWPOS-0110-2)
    (new-gname 'WWPOS-0110-3)
    (new-gname 'WWPOS-0101-2)
    (new-gname 'WWPOS-0101-4)
    (new-gname 'WWPOS-0011-3)
    (new-gname 'WWPOS-0011-4)
    (new-gname 'WWPOS-0111-2)
    (new-gname 'WWPOS-0111-3)
    (new-gname 'WWPOS-0111-4)
    (new-gname 'WWPOS-1011-1)
    (new-gname 'WWPOS-1011-3)
    (new-gname 'WWPOS-1011-4)
    (new-gname 'WWPOS-1101-1)
    (new-gname 'WWPOS-1101-2)
    (new-gname 'WWPOS-1101-4)
    (new-gname 'WWPOS-1110-1)
    (new-gname 'WWPOS-1110-2)
    (new-gname 'WWPOS-1110-3)
    (new-gname 'WWPOS-1111-1)
    (new-gname 'WWPOS-1111-2)
    (new-gname 'WWPOS-1111-3)
    (new-gname 'WWPOS-1111-4)))
@

The order of these named lists is critical and cannot be changed!
For any state in new-fast-latent-states etc, the corresponding element of 
   this array records what strain they have (it must correspond to the
   final digit of the symbol name).

<<>>=
(defvar major-strain
  (list 1 2 3 4 1 2 1 3 1 4 2 3 2 4 3 4 2 3 4 1 3 4 1 2 4 1 2 3 1 2 3 4))

; Added SA 092112
(defvar inh-sensitive-dom-p
  (list 1 0 1 0 1 0 1 1 1 0 0 1 0 0 1 0 0 1 0 1 1 0 1 0 0 1 0 1 1 0 1 0))

(defvar inh-sensitive-p    ; if you take isoniazid will it kill everything you have
  (list 1 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0))

; Added SA 092012
; Dominant strain not resistant to rifampicin => fails standard treatment
(defvar rif-sensitive-n   
  (list 0 0 1 1 0 0 0 1 0 1 0 1 0 1 1 1 0 1 1 0 1 1 0 0 1 0 0 1 0 0 1 1))
; End edits SA 092012
@

<<>>=
(defvar tb-state-list 
  (append
   (list (new-gname 'UU) (new-gname 'ZZ) (new-gname 'QQ))
   new-fast-latent-states  
   new-slow-latent-states  
   smearnegative-states
   smearpositive-states
   old-fast-latent-states
   old-slow-latent-states
   incomplete-cure-states
   treated-smneg-states-standard
   treated-smneg-states-specialized
   treated-smpos-states-standard
   treated-smpos-states-specialized))
@

<<>>=
(add-state-list tbmodel tb-state-list)

(add-state tbmodel (new-gname 'cumul-tb-mort-before-treatment))
(add-state tbmodel (new-gname 'cumul-tb-mort-after-treatment))
@

Progression of fast latent infections:
  Assumption: the rate of progression can depend on drug resistance categories
  However: the probability of smear positivity does not.
<<>>=
(defvar gammas (list (new-gname 'gamma1) (new-gname 'gamma2) (new-gname 'gamma3) (new-gname 'gamma4)))
(add-parameter tbmodel (new-gname 'gamma1))
(add-parameter tbmodel (new-gname 'gamma2))
(add-parameter tbmodel (new-gname 'gamma3))
(add-parameter tbmodel (new-gname 'gamma4))
(add-parameter tbmodel (new-gname 'smpos-prob))
(defvar ns (length new-fast-latent-states))
(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii new-fast-latent-states)
                           :to (nth ii smearnegative-states)
                           :rate (s* (nth ii new-fast-latent-states)
                                     (nth (- (nth ii major-strain) 1) gammas)
                                     (s- 1 (new-gname 'smpos-prob)))) 
  (add-transition tbmodel :from (nth ii new-fast-latent-states)
                           :to (nth ii smearpositive-states)
                           :rate (s* (nth ii new-fast-latent-states)
                                     (nth (- (nth ii major-strain) 1) gammas)
                                     (new-gname 'smpos-prob))))
(setf tbmodel1 (deep-copy-model tbmodel))
@  

Progression of slow latent infections:
  Assumption: the rate of progression can depend on drug resistance categories
  However: the probability of smear positivity does not.
Also note: the smear positive probability isn't depending on fast/slow either
<<>>=
(defvar gammaprimes 
  (list (new-gname 'gamma1prime) (new-gname 'gamma2prime) (new-gname 'gamma3prime) (new-gname 'gamma4prime)))
(add-parameter tbmodel (new-gname 'gamma1prime))
(add-parameter tbmodel (new-gname 'gamma2prime))
(add-parameter tbmodel (new-gname 'gamma3prime))
(add-parameter tbmodel (new-gname 'gamma4prime))
(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii new-slow-latent-states)
                          :to (nth ii smearnegative-states)
                          :rate (s* (nth ii new-slow-latent-states)
                                    (nth (- (nth ii major-strain) 1) gammaprimes)
                                    (s- 1 (new-gname 'smpos-prob)))) 
  (add-transition tbmodel :from (nth ii new-slow-latent-states)
                          :to (nth ii smearpositive-states)
                          :rate (s* (nth ii new-slow-latent-states)
                                    (nth (- (nth ii major-strain) 1) gammaprimes)
                                    (new-gname 'smpos-prob))))
@

Fast latency to slow latency:
Assumption: this does not depend on any drug resistant state.
<<>>=
(add-parameter tbmodel (new-gname 'eta))
(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii new-fast-latent-states)
                          :to (nth ii new-slow-latent-states)
                          :rate (s* (nth ii new-fast-latent-states)
                                    (new-gname 'eta))))

(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii old-fast-latent-states)
                           :to (nth ii smearnegative-states)
                           :rate (s* (nth ii old-fast-latent-states)
                                     (nth (- (nth ii major-strain) 1) gammas)
                                     (s- 1 (new-gname 'smpos-prob)))) 
  (add-transition tbmodel :from (nth ii old-fast-latent-states)
                           :to (nth ii smearpositive-states)
                           :rate (s* (nth ii old-fast-latent-states)
                                     (nth (- (nth ii major-strain) 1) gammas)
                                     (new-gname 'smpos-prob))))
  
(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii old-slow-latent-states)
                          :to (nth ii smearnegative-states)
                          :rate (s* (nth ii old-slow-latent-states)
                                    (nth (- (nth ii major-strain) 1) gammaprimes)
                                    (s- 1 (new-gname 'smpos-prob)))) 
  (add-transition tbmodel :from (nth ii old-slow-latent-states)
                          :to (nth ii smearpositive-states)
                          :rate (s* (nth ii old-slow-latent-states)
                                    (nth (- (nth ii major-strain) 1) gammaprimes)
                                    (new-gname 'smpos-prob))))
@

Fast latency to slow latency.
Assumption: this does not depend on any drug resistant state.
<<>>=
(add-parameter tbmodel (new-gname 'eta))
(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii old-fast-latent-states)
                          :to (nth ii old-slow-latent-states)
                          :rate (s* (nth ii old-fast-latent-states)
                                    (new-gname 'eta))))
@

Chemoprophylaxis/LTBI therapy
<<>>=
(add-parameter tbmodel (new-gname 'sigma))
(add-parameter tbmodel (new-gname 'sigma-old))    ; chemo rate for ltbi in ex-cases
(dotimes (ii ns)
  (if (equal (nth ii inh-sensitive-p) 1)
      (add-transition tbmodel :from (nth ii new-slow-latent-states)
                              :to (new-gname 'qq)
                              :rate (s* (nth ii new-slow-latent-states)
                                        (new-gname 'sigma)))))
(dotimes (ii ns)
  (if (equal (nth ii inh-sensitive-p) 1)
      (add-transition tbmodel :from (nth ii new-fast-latent-states)
                              :to (new-gname 'qq)
                              :rate (s* (nth ii new-fast-latent-states)
                                        (new-gname 'sigma)))))


(dotimes (ii ns)
  (if (equal (nth ii inh-sensitive-p) 1)
      (add-transition tbmodel :from (nth ii old-slow-latent-states)
                              :to (new-gname 'zz)
                              :rate (s* (nth ii old-slow-latent-states)
                                        (new-gname 'sigma-old)))))
(dotimes (ii ns)
  (if (equal (nth ii inh-sensitive-p) 1)
      (add-transition tbmodel :from (nth ii old-fast-latent-states)
                              :to (new-gname 'zz)
                              :rate (s* (nth ii old-fast-latent-states)
                                        (new-gname 'sigma-old)))))
@

Smear positive to negative, undiagnosed:
<<>>=
(add-parameter tbmodel (new-gname 'theta))
(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii smearpositive-states)
                          :to (nth ii smearnegative-states)
                          :rate (s* (nth ii smearpositive-states)
                                    (new-gname 'theta))))
@

Smear negative to positive, undiagnosed:
<<>>=
(add-parameter tbmodel (new-gname 'kappa))
(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii smearnegative-states)
                          :to (nth ii smearpositive-states)
                          :rate (s* (nth ii smearnegative-states)
                                    (new-gname 'kappa))))
@

Self heal:
does NOT depend on drug resistance.
<<>>=
(add-parameter tbmodel (new-gname 'rhoprime))
(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii smearnegative-states)
                          :to (nth ii incomplete-cure-states)
                          :rate (s* (nth ii smearnegative-states)
                                    (new-gname 'rhoprime))))
(add-parameter tbmodel (new-gname 'rho))
(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii smearpositive-states)
                          :to (nth ii incomplete-cure-states)
                          :rate (s* (nth ii smearpositive-states)
                                    (new-gname 'rho))))
@

Relapse after self-heal (without treatment):
<<>>=
; Added SA 092012
(defvar phis (list (new-gname 'phi-wild) (new-gname 'phi-inh) (new-gname 'phi-rif) (new-gname 'phi-mdr)))
(add-parameter tbmodel (new-gname 'phi-wild))
(add-parameter tbmodel (new-gname 'phi-inh))
(add-parameter tbmodel (new-gname 'phi-rif))
(add-parameter tbmodel (new-gname 'phi-mdr))
(add-parameter tbmodel (new-gname 'smpos-prob))
(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii incomplete-cure-states)
                          :to (nth ii smearnegative-states)
                          :rate (s* (nth ii incomplete-cure-states)
                                    (nth (- (nth ii major-strain) 1) phis)
				    (s- 1 (new-gname 'smpos-prob)))))

(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii incomplete-cure-states)
                          :to (nth ii smearpositive-states)
                          :rate (s* (nth ii incomplete-cure-states)
                                    (nth (- (nth ii major-strain) 1) phis)
				    (new-gname 'smpos-prob))))
; End SA edit 092012
@

<<>>=
;; SA edit 101512
(add-parameter tbmodel (new-gname 'omega))
(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii incomplete-cure-states)
                          :to (nth ii old-slow-latent-states)
                          :rate (s* (nth ii incomplete-cure-states)
                                    (new-gname 'omega))))
;; End SA edit 101512
@

Treatment
<<>>=
;   diagnosis rate could depend on smear status
; SA 092012 -- Changed delta to py for consistency with latex eqns
(add-parameter tbmodel (new-gname 'pyprime))
(add-parameter tbmodel (new-gname 'drugtest1))
(add-parameter tbmodel (new-gname 'drugtest2))
(add-parameter tbmodel (new-gname 'drugtest3))
(add-parameter tbmodel (new-gname 'drugtest4))
(setf drugtest-results (list (new-gname 'drugtest1) (new-gname 'drugtest2) (new-gname 'drugtest3) (new-gname 'drugtest4)))
(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii smearnegative-states)
                          :to (nth ii treated-smneg-states-specialized)
                          :rate (s* (nth ii smearnegative-states)
                                    (new-gname 'pyprime)
                                    (nth (- (nth ii major-strain) 1) drugtest-results)))
  (add-transition tbmodel :from (nth ii smearnegative-states)
                          :to (nth ii treated-smneg-states-standard)
                          :rate (s* (nth ii smearnegative-states)
                                    (new-gname 'pyprime) 
                                    (s- 1 (nth (- (nth ii major-strain) 1) drugtest-results)))))
@

<<>>=
; SA 092012 -- Changed delta to py for consistency with latex eqns
(add-parameter tbmodel (new-gname 'py))
(dotimes (ii ns)
  (setf tmp-majorstrain (nth ii major-strain))
  (setf inh-p (nth ii inh-sensitive-dom-p))
  (if inh-p (setf tmp-finalstates (list (nth ii treated-smneg-states-standard) (nth ii treated-smneg-states-specialized)))
    (setf tmp-finalstates (list (nth ii treated-smpos-states-standard) (nth ii treated-smpos-states-specialized))))
  (add-transition tbmodel :from (nth ii smearpositive-states)
                          :to (cadr tmp-finalstates)
                          :rate (s* (nth ii smearpositive-states)
                                    (new-gname 'py)
                                    (nth (- tmp-majorstrain 1) drugtest-results)))
  (add-transition tbmodel :from (nth ii smearpositive-states)
                          :to (car tmp-finalstates)
                          :rate (s* (nth ii smearpositive-states)
                                    (new-gname 'py) 
                                    (s- 1 (nth (- tmp-majorstrain 1) drugtest-results)))))

<<>>=
; Added 092012 SA 
; @@@
; For some treatment ends without cure--return to being active cases
(add-parameter tbmodel (new-gname 'treatment-end-no-cure-yy))
(dotimes (ii ns)
  (cond ((equal (nth ii rif-sensitive-n) 1)
         (add-transition tbmodel :from (nth ii treated-smpos-states-standard)
	                          :to (nth ii smearpositive-states)
			          
			          :rate (s* (nth ii treated-smpos-states-standard)
				            (new-gname 'treatment-end-no-cure-yy)))
          (add-transition tbmodel :from (nth ii treated-smneg-states-standard)
                                  :to (nth ii smearnegative-states)
     			          
     			          :rate (s* (nth ii treated-smneg-states-standard)
				            (new-gname 'treatment-end-no-cure-yy))))))
; End edits SA 092012
@

Smear positives become smearnegatives while on treatment.
<<>>=
(defvar taus (list (new-gname 'tau1) (new-gname 'tau2) (new-gname 'tau3) (new-gname 'tau4)))
(add-parameter tbmodel (new-gname 'tau1))
(add-parameter tbmodel (new-gname 'tau2))
(add-parameter tbmodel (new-gname 'tau3))
(add-parameter tbmodel (new-gname 'tau4))

(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii treated-smpos-states-standard)
                          :to (nth ii treated-smneg-states-standard)
                          :rate (s* (nth ii treated-smpos-states-standard)
                                    (nth (- (nth ii major-strain) 1) taus))))

(defvar tauprimes (list (new-gname 'tau1prime) (new-gname 'tau2prime) (new-gname 'tau3prime) (new-gname 'tau4prime)))
(add-parameter tbmodel (new-gname 'tau1prime))
(add-parameter tbmodel (new-gname 'tau2prime))
(add-parameter tbmodel (new-gname 'tau3prime))
(add-parameter tbmodel (new-gname 'tau4prime))
(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii treated-smpos-states-specialized)
                          :to (nth ii treated-smneg-states-specialized)
                          :rate (s* (nth ii treated-smpos-states-specialized)
                                    (nth (- (nth ii major-strain) 1) tauprimes))))  
@

During treatment, some incorrectly on standard treatment move to specialized treatment:
<<>>=
(defvar nupos (list (new-gname 'nu1pos) (new-gname 'nu2pos) (new-gname 'nu3pos) (new-gname 'nu4pos)))
(add-parameter tbmodel (new-gname 'nu1pos))
(add-parameter tbmodel (new-gname 'nu2pos))
(add-parameter tbmodel (new-gname 'nu3pos))
(add-parameter tbmodel (new-gname 'nu4pos))
@

<<>>=
(defvar nuneg (list (new-gname 'nu1neg) (new-gname 'nu2neg) (new-gname 'nu3neg) (new-gname 'nu4neg)))
(add-parameter tbmodel (new-gname 'nu1neg))
(add-parameter tbmodel (new-gname 'nu2neg))
(add-parameter tbmodel (new-gname 'nu3neg))
(add-parameter tbmodel (new-gname 'nu4neg))
(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii treated-smpos-states-standard)
                          :to (nth ii treated-smpos-states-specialized)
                          :rate (s* (nth ii treated-smpos-states-standard)
                                    (nth (- (nth ii major-strain) 1) nupos))))
(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii treated-smneg-states-standard)
                          :to (nth ii treated-smneg-states-specialized)
                          :rate (s* (nth ii treated-smneg-states-standard)
                                    (nth (- (nth ii major-strain) 1) nuneg))))

; @@@
;(setf tbmodel1 (deep-copy-model tbmodel))
; People fail treatment 

<<>>=
(add-parameter tbmodel (new-gname 'treatment-interrupt))
(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii treated-smneg-states-standard)
                          :to (nth ii smearnegative-states)
                          :rate (s* (nth ii treated-smneg-states-standard)
                                    (new-gname 'treatment-interrupt))))

(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii treated-smneg-states-specialized)
                          :to (nth ii smearnegative-states)
                          :rate (s* (nth ii treated-smneg-states-specialized)
                                    (new-gname 'treatment-interrupt))))

(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii treated-smpos-states-standard)
                          :to (nth ii smearpositive-states)
                          :rate (s* (nth ii treated-smpos-states-standard)
                                    (new-gname 'treatment-interrupt))))

(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii treated-smpos-states-specialized)
                          :to (nth ii smearpositive-states)
                          :rate (s* (nth ii treated-smpos-states-specialized)
                                    (new-gname 'treatment-interrupt))))
@

TB death:
<<>>=
(add-parameter tbmodel (new-gname 'smpos-mort-untrt))
(dotimes (ii ns) 
  (add-transition tbmodel :from (nth ii smearpositive-states)
                          :to (new-gname 'cumul-tb-mort-before-treatment)
                          :rate (s* (nth ii smearpositive-states) (new-gname 'smpos-mort-untrt))))
(add-parameter tbmodel (new-gname 'smneg-mort-untrt))
(dotimes (ii ns) 
  (add-transition tbmodel :from (nth ii smearnegative-states)
                          :to (new-gname 'cumul-tb-mort-before-treatment)
                          :rate (s* (nth ii smearnegative-states) (new-gname 'smneg-mort-untrt))))
@

<<>>=
; right now mortality depends on nothing
;; Edits SA 101512
(defvar tb-mort-trt-yy (list (new-gname 'mort-treat-wild-yy) (new-gname 'mort-treat-inh-yy) (new-gname 'mort-treat-rif-yy) (new-gname 'mort-treat-mdr-yy)))

(add-parameter tbmodel (new-gname 'mort-treat-wild-yy))
(add-parameter tbmodel (new-gname 'mort-treat-inh-yy))
(add-parameter tbmodel (new-gname 'mort-treat-rif-yy))
(add-parameter tbmodel (new-gname 'mort-treat-mdr-yy))

(defvar tb-mort-trt-ww (list (new-gname 'mort-treat-wild-ww) (new-gname 'mort-treat-inh-ww) (new-gname 'mort-treat-rif-ww) (new-gname 'mort-treat-mdr-ww)))

(add-parameter tbmodel (new-gname 'mort-treat-wild-ww))
(add-parameter tbmodel (new-gname 'mort-treat-inh-ww))
(add-parameter tbmodel (new-gname 'mort-treat-rif-ww))
(add-parameter tbmodel (new-gname 'mort-treat-mdr-ww))

(dotimes (ii ns) 
  (add-transition tbmodel :from (nth ii treated-smneg-states-standard)
                          :to (new-gname 'cumul-tb-mort-after-treatment)
                          :rate (s* (nth ii treated-smneg-states-standard) 
				    (nth (- (nth ii major-strain) 1) tb-mort-trt-yy))))

(dotimes (ii ns) 
  (add-transition tbmodel :from (nth ii treated-smneg-states-specialized)
                          :to (new-gname 'cumul-tb-mort-after-treatment)
                          :rate (s* (nth ii treated-smneg-states-specialized) 
				    (nth (- (nth ii major-strain) 1) tb-mort-trt-ww))))

(dotimes (ii ns) 
  (add-transition tbmodel :from (nth ii treated-smpos-states-standard)
                          :to (new-gname 'cumul-tb-mort-after-treatment)
                          :rate (s* (nth ii treated-smpos-states-standard) 
				    (nth (- (nth ii major-strain) 1) tb-mort-trt-yy))))

(dotimes (ii ns) 
  (add-transition tbmodel :from (nth ii treated-smpos-states-specialized)
                          :to (new-gname 'cumul-tb-mort-after-treatment)
                          :rate (s* (nth ii treated-smpos-states-specialized) 
				    (nth (- (nth ii major-strain) 1) tb-mort-trt-ww))))

;; End edits SA 101215
@

Transitions related to resistance and treatment:
The key phenomena here are 
\begin{itemize}
\item the sequential addition of resistance during treatment, with our assumption being that 
    it is easier to add further resistance once the patient is already acquired some drug resistance
\item For now, we do not distinguish between acquiring INH resistance and RIF resistance, but we might
\item Nearly all patients who are treated (who do not die) achieve complete cure and are no longer
    at risk for tuberculosis, but a small fraction remain at risk of relapse (RR)
\item Cure takes a different amount of time depending on the level of drug resistance, with 
    cure-wild-yy (time) < cure-inh-yy << cure-rif-yy < cure-mdr-yy
    We are not distinguishing XDR...yet.
\end{itemize}
 
 Our assumptions are that the predominant strain has enough population diversity to be able to
   select mutations for further resistance.  New strains can endogenously arise from the dominant
   strain only, and only sequentially.  1WT -> 2INHR or 3RIFR.  2INHR-> 4MDR.  3RIFR->4MDR.
 During treatment, resistant strains can overgrow less resistant strains, and we have not
   distinguished these probabilities

Subscripting changing transitions:
People are cured, people develop drug resistance, standard treatment 

<<>>=
;; Edited 101912 SA
;;

(defvar comp-yys (list (new-gname 'comp-rifsen-yy) (new-gname 'comp-rifres-yy)))
(add-parameter tbmodel (new-gname 'comp-rifsen-yy))
(add-parameter tbmodel (new-gname 'comp-rifres-yy))

(defvar cure-yys (list (new-gname 'cure-wild-yy) (new-gname 'cure-inh-yy) (new-gname 'cure-rif-yy) (new-gname 'cure-mdr-yy)))
(add-parameter tbmodel (new-gname 'cure-wild-yy))
(add-parameter tbmodel (new-gname 'cure-inh-yy))
(add-parameter tbmodel (new-gname 'cure-rif-yy))
(add-parameter tbmodel (new-gname 'cure-mdr-yy))

(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii treated-smneg-states-standard)
		  :to (new-gname 'ZZ)                          
		  
		  :rate (s* (nth ii treated-smneg-states-standard)
			    (nth (nth ii rif-sensitive-n) comp-yys)
			    (nth (- (nth ii major-strain) 1) cure-yys))))

(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii treated-smneg-states-standard)
		  :to (nth ii incomplete-cure-states)                          
		  
		  :rate (s* (nth ii treated-smneg-states-standard)
			    (s- 1 (nth (nth ii rif-sensitive-n) comp-yys))
			    (nth (- (nth ii major-strain) 1) cure-yys))))
@

<<>>=
(defvar comp-wws (list (new-gname 'comp-rifsen-ww) (new-gname 'comp-rifres-ww)))
(add-parameter tbmodel (new-gname 'comp-rifsen-ww))
(add-parameter tbmodel (new-gname 'comp-rifres-ww))

(defvar cure-wws (list (new-gname 'cure-wild-ww) (new-gname 'cure-inh-ww) (new-gname 'cure-rif-ww) (new-gname 'cure-mdr-ww)))
(add-parameter tbmodel (new-gname 'cure-wild-ww))
(add-parameter tbmodel (new-gname 'cure-inh-ww))
(add-parameter tbmodel (new-gname 'cure-rif-ww))
(add-parameter tbmodel (new-gname 'cure-mdr-ww))

(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii treated-smneg-states-specialized)
		  :to (new-gname 'ZZ)                          
		  
		  :rate (s* (nth ii treated-smneg-states-standard)
			    (nth (nth ii rif-sensitive-n) comp-wws)
			    (nth (- (nth ii major-strain) 1) cure-wws))))

(dotimes (ii ns)
  (add-transition tbmodel :from (nth ii treated-smneg-states-specialized)
		  :to (nth ii incomplete-cure-states)                          
		  
		  :rate (s* (nth ii treated-smneg-states-standard)
			    (s- 1 (nth (nth ii rif-sensitive-n) comp-wws))
			    (nth (- (nth ii major-strain) 1) cure-wws))))
@

<<>>=
(add-parameter tbmodel (new-gname 'strainacc-inh-yy))
(add-parameter tbmodel (new-gname 'strainacc-rif-yy))
(add-parameter tbmodel (new-gname 'strainacc-mdr-yy))
(add-parameter tbmodel (new-gname 'overgrowth-yy))
;; End SA edits 101512/101912
@

People are cured, people develop drug resistance, specialized treatment
<<>>=
(add-parameter tbmodel (new-gname 'strainacc-rif-ww))
(add-parameter tbmodel (new-gname 'strainacc-inh-ww))
(add-parameter tbmodel (new-gname 'strainacc-mdr-ww))
(add-parameter tbmodel (new-gname 'overgrowth-ww))

; WWPOS

(add-transition tbmodel  
		 
		:from  (new-gname 'WWPOS-1000-1) 
		:to (new-gname 'WWPOS-1100-1) 
		:rate (s* (new-gname 'WWPOS-1000-1) (new-gname 'strainacc-inh-ww)))

(add-transition tbmodel  
		 
		:from  (new-gname 'WWPOS-1000-1) 
		:to (new-gname 'WWPOS-1010-1) 
		:rate (s* (new-gname 'WWPOS-1000-1) (new-gname 'strainacc-rif-ww)))

(add-transition tbmodel  
		 
		:from (new-gname 'WWPOS-0100-2) 
		:to (new-gname 'WWPOS-0101-2) 
		:rate (s* (new-gname 'WWPOS-0100-2) (new-gname 'strainacc-mdr-ww)))


(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-0010-3) 
	:to (new-gname 'WWPOS-0011-3) 
	:rate (s* (new-gname 'WWPOS-0010-3) (new-gname 'strainacc-mdr-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1100-1) 
	:to (new-gname 'WWPOS-1110-1) 
	:rate (s* (new-gname 'WWPOS-1100-1) (new-gname 'strainacc-rif-ww)))
(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1100-1) 
	:to (new-gname 'WWPOS-1100-2) 
	:rate (s* (new-gname 'WWPOS-1100-1) (new-gname 'overgrowth-ww)))


(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1100-2) 
	:to (new-gname 'WWPOS-1101-2) 
	:rate (s* (new-gname 'WWPOS-1100-2) (new-gname 'strainacc-mdr-ww)))
@
<<>>=
(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1010-1) 
	:to (new-gname 'WWPOS-1110-1) 
	:rate (s* (new-gname 'WWPOS-1010-1) (new-gname 'strainacc-inh-ww)))
(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1010-1) 
	:to (new-gname 'WWPOS-1010-3) 
	:rate (s* (new-gname 'WWPOS-1010-1) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1010-3) 
	:to (new-gname 'WWPOS-1011-3) 
	:rate (s* (new-gname 'WWPOS-1010-3) (new-gname 'strainacc-mdr-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1001-1) 
	:to (new-gname 'WWPOS-1101-1) 
	:rate (s* (new-gname 'WWPOS-1001-1) (new-gname 'strainacc-inh-ww)))
(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1001-1) 
	:to (new-gname 'WWPOS-1011-1) 
	:rate (s* (new-gname 'WWPOS-1001-1) (new-gname 'strainacc-rif-ww)))
(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1001-1) 
	:to (new-gname 'WWPOS-1001-4) 
	:rate (s* (new-gname 'WWPOS-1001-1) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-0110-2) 
	:to (new-gname 'WWPOS-0111-2) 
	:rate (s* (new-gname 'WWPOS-0110-2) (new-gname 'strainacc-mdr-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-0110-3) 
	:to (new-gname 'WWPOS-0111-3) 
	:rate (s* (new-gname 'WWPOS-0110-3) (new-gname 'strainacc-mdr-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-0101-2) 
	:to (new-gname 'WWPOS-0101-4) 
	:rate (s* (new-gname 'WWPOS-0101-2) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-0011-3) 
	:to (new-gname 'WWPOS-0011-4) 
	:rate (s* (new-gname 'WWPOS-0011-3) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1110-1) 
	:to (new-gname 'WWPOS-1110-2) 
	:rate (s* (new-gname 'WWPOS-1110-1) (new-gname 'overgrowth-ww) ))
(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1110-1) 
	:to (new-gname 'WWPOS-1110-3) 
	:rate (s* (new-gname 'WWPOS-1110-1) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1110-2) 
	:to (new-gname 'WWPOS-1111-2) 
	:rate (s* (new-gname 'WWPOS-1110-2) (new-gname 'strainacc-mdr-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1110-3) 
	:to (new-gname 'WWPOS-1111-3) 
	:rate (s* (new-gname 'WWPOS-1110-3) (new-gname 'strainacc-mdr-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1101-1) 
	:to (new-gname 'WWPOS-1111-1) 
	:rate (s* (new-gname 'WWPOS-1101-1) (new-gname 'strainacc-rif-ww) ))
(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1101-1) 
	:to (new-gname 'WWPOS-1101-2) 
	:rate (s* (new-gname 'WWPOS-1101-1) (new-gname 'overgrowth-ww)))
(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1101-1) 
	:to (new-gname 'WWPOS-1101-4) 
	:rate (s* (new-gname 'WWPOS-1101-1) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1101-2) 
	:to (new-gname 'WWPOS-1101-4) 
	:rate (s* (new-gname 'WWPOS-1101-2) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1011-1) 
	:to (new-gname 'WWPOS-1111-1) 
	:rate (s* (new-gname 'WWPOS-1011-1) (new-gname 'strainacc-inh-ww)))
(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1011-1) 
	:to (new-gname 'WWPOS-1011-3) 
	:rate (s* (new-gname 'WWPOS-1011-1) (new-gname 'overgrowth-ww)))
(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1011-1) 
	:to (new-gname 'WWPOS-1011-4) 
	:rate (s* (new-gname 'WWPOS-1011-1) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1011-3) 
	:to (new-gname 'WWPOS-1011-4) 
	:rate (s* (new-gname 'WWPOS-1011-3) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-0111-2) 
	:to (new-gname 'WWPOS-0111-4) 
	:rate (s* (new-gname 'WWPOS-0111-2) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-0111-3) 
	:to (new-gname 'WWPOS-0111-4) 
	:rate (s* (new-gname 'WWPOS-0111-3) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1111-1) 
	:to (new-gname 'WWPOS-1111-2) 
	:rate (s* (new-gname 'WWPOS-1111-1) (new-gname 'overgrowth-ww)))
(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1111-1) 
	:to (new-gname 'WWPOS-1111-3) 
	:rate (s* (new-gname 'WWPOS-1111-1) (new-gname 'overgrowth-ww)))
(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1111-1) 
	:to (new-gname 'WWPOS-1111-4) 
	:rate (s* (new-gname 'WWPOS-1111-1) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1111-2) 
	:to (new-gname 'WWPOS-1111-4) 
	:rate (s* (new-gname 'WWPOS-1111-2) (new-gname 'overgrowth-ww) ))

(add-transition tbmodel  
	 
	:from (new-gname 'WWPOS-1111-3) 
	:to (new-gname 'WWPOS-1111-4) 
	:rate (s* (new-gname 'WWPOS-1111-3) (new-gname 'overgrowth-ww) ))
@

<<>>=
(add-transition tbmodel  
	 
	:from  (new-gname 'WWNEG-1000-1) 
	:to (new-gname 'WWNEG-1100-1) 
	:rate (s* (new-gname 'WWNEG-1000-1) (new-gname 'strainacc-inh-ww)))
(add-transition tbmodel  
	 
	:from  (new-gname 'WWNEG-1000-1) 
	:to (new-gname 'WWNEG-1010-1) 
	:rate (s* (new-gname 'WWNEG-1000-1) (new-gname 'strainacc-rif-ww)))


(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-0100-2) 
	:to (new-gname 'WWNEG-0101-2) 
	:rate (s* (new-gname 'WWNEG-0100-2) (new-gname 'strainacc-mdr-ww)))


(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-0010-3) 
	:to (new-gname 'WWNEG-0011-3) 
	:rate (s* (new-gname 'WWNEG-0010-3) (new-gname 'strainacc-mdr-ww)))


(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1100-1) 
	:to (new-gname 'WWNEG-1110-1) 
	:rate (s* (new-gname 'WWNEG-1100-1) (new-gname 'strainacc-rif-ww)))
(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1100-1) 
	:to (new-gname 'WWNEG-1100-2) 
	:rate (s* (new-gname 'WWNEG-1100-1) (new-gname 'overgrowth-ww)))


(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1100-2) 
	:to (new-gname 'WWNEG-1101-2) 
	:rate (s* (new-gname 'WWNEG-1100-2) (new-gname 'strainacc-mdr-ww)))


(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1010-1) 
	:to (new-gname 'WWNEG-1110-1) 
	:rate (s* (new-gname 'WWNEG-1010-1) (new-gname 'strainacc-inh-ww)))
(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1010-1) 
	:to (new-gname 'WWNEG-1010-3) 
	:rate (s* (new-gname 'WWNEG-1010-1) (new-gname 'overgrowth-ww)))


(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1010-3) 
	:to (new-gname 'WWNEG-1011-3) 
	:rate (s* (new-gname 'WWNEG-1010-3) (new-gname 'strainacc-mdr-ww)))


(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1001-1) 
	:to (new-gname 'WWNEG-1101-1) 
	:rate (s* (new-gname 'WWNEG-1001-1) (new-gname 'strainacc-inh-ww)))
(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1001-1) 
	:to (new-gname 'WWNEG-1011-1) 
	:rate (s* (new-gname 'WWNEG-1001-1) (new-gname 'strainacc-rif-ww)))
(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1001-1) 
	:to (new-gname 'WWNEG-1001-4) 
	:rate (s* (new-gname 'WWNEG-1001-1) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-0110-2) 
	:to (new-gname 'WWNEG-0111-2) 
	:rate (s* (new-gname 'WWNEG-0110-2) (new-gname 'strainacc-mdr-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-0110-3) 
	:to (new-gname 'WWNEG-0111-3) 
	:rate (s* (new-gname 'WWNEG-0110-3) (new-gname 'strainacc-mdr-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-0101-2) 
	:to (new-gname 'WWNEG-0101-4) 
	:rate (s* (new-gname 'WWNEG-0101-2) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-0011-3) 
	:to (new-gname 'WWNEG-0011-4) 
	:rate (s* (new-gname 'WWNEG-0011-3) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1110-1) 
	:to (new-gname 'WWNEG-1110-2) 
	:rate (s* (new-gname 'WWNEG-1110-1) (new-gname 'overgrowth-ww) ))
(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1110-1) 
	:to (new-gname 'WWNEG-1110-3) 
	:rate (s* (new-gname 'WWNEG-1110-1) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1110-2) 
	:to (new-gname 'WWNEG-1111-2) 
	:rate (s* (new-gname 'WWNEG-1110-2) (new-gname 'strainacc-mdr-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1110-3) 
	:to (new-gname 'WWNEG-1111-3) 
	:rate (s* (new-gname 'WWNEG-1110-3) (new-gname 'strainacc-mdr-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1101-1) 
	:to (new-gname 'WWNEG-1111-1) 
	:rate (s* (new-gname 'WWNEG-1101-1) (new-gname 'strainacc-rif-ww) ))
(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1101-1) 
	:to (new-gname 'WWNEG-1101-2) 
	:rate (s* (new-gname 'WWNEG-1101-1) (new-gname 'overgrowth-ww)))
(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1101-1) 
	:to (new-gname 'WWNEG-1101-4) 
	:rate (s* (new-gname 'WWNEG-1101-1) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1101-2) 
	:to (new-gname 'WWNEG-1101-4) 
	:rate (s* (new-gname 'WWNEG-1101-2) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1011-1) 
	:to (new-gname 'WWNEG-1111-1) 
	:rate (s* (new-gname 'WWNEG-1011-1) (new-gname 'strainacc-inh-ww)))
(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1011-1) 
	:to (new-gname 'WWNEG-1011-3) 
	:rate (s* (new-gname 'WWNEG-1011-1) (new-gname 'overgrowth-ww)))
(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1011-1) 
	:to (new-gname 'WWNEG-1011-4) 
	:rate (s* (new-gname 'WWNEG-1011-1) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1011-3) 
	:to (new-gname 'WWNEG-1011-4) 
	:rate (s* (new-gname 'WWNEG-1011-3) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-0111-2) 
	:to (new-gname 'WWNEG-0111-4) 
	:rate (s* (new-gname 'WWNEG-0111-2) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-0111-3) 
	:to (new-gname 'WWNEG-0111-4) 
	:rate (s* (new-gname 'WWNEG-0111-3) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1111-1) 
	:to (new-gname 'WWNEG-1111-2) 
	:rate (s* (new-gname 'WWNEG-1111-1) (new-gname 'overgrowth-ww)))
(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1111-1) 
	:to (new-gname 'WWNEG-1111-3) 
	:rate (s* (new-gname 'WWNEG-1111-1) (new-gname 'overgrowth-ww)))
(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1111-1) 
	:to (new-gname 'WWNEG-1111-4) 
	:rate (s* (new-gname 'WWNEG-1111-1) (new-gname 'overgrowth-ww)))

(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1111-2) 
	:to (new-gname 'WWNEG-1111-4) 
	:rate (s* (new-gname 'WWNEG-1111-2) (new-gname 'overgrowth-ww) ))


(add-transition tbmodel  
	 
	:from (new-gname 'WWNEG-1111-3) 
	:to (new-gname 'WWNEG-1111-4) 
	:rate (s* (new-gname 'WWNEG-1111-3) (new-gname 'overgrowth-ww) ))
@

YYPOS
<<>>=
(add-transition tbmodel  
	 
	:from  (new-gname 'YYPOS-1000-1) 
	:to (new-gname 'YYPOS-1100-1) 
	:rate (s* (new-gname 'YYPOS-1000-1) (new-gname 'strainacc-inh-yy)))
(add-transition tbmodel  
	 
	:from  (new-gname 'YYPOS-1000-1) 
	:to (new-gname 'YYPOS-1010-1) 
	:rate (s* (new-gname 'YYPOS-1000-1) (new-gname 'strainacc-rif-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-0100-2) 
	:to (new-gname 'YYPOS-0101-2) 
	:rate (s* (new-gname 'YYPOS-0100-2) (new-gname 'strainacc-mdr-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-0010-3) 
	:to (new-gname 'YYPOS-0011-3) 
	:rate (s* (new-gname 'YYPOS-0010-3) (new-gname 'strainacc-mdr-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1100-1) 
	:to (new-gname 'YYPOS-1110-1) 
	:rate (s* (new-gname 'YYPOS-1100-1) (new-gname 'strainacc-rif-yy)))
(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1100-1) 
	:to (new-gname 'YYPOS-1100-2) 
	:rate (s* (new-gname 'YYPOS-1100-1) (new-gname 'overgrowth-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1100-2) 
	:to (new-gname 'YYPOS-1101-2) 
	:rate (s* (new-gname 'YYPOS-1100-2) (new-gname 'strainacc-mdr-yy)))


(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1010-1) 
	:to (new-gname 'YYPOS-1110-1) 
	:rate (s* (new-gname 'YYPOS-1010-1) (new-gname 'strainacc-inh-yy)))
(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1010-1) 
	:to (new-gname 'YYPOS-1010-3) 
	:rate (s* (new-gname 'YYPOS-1010-1) (new-gname 'overgrowth-yy)))


(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1010-3) 
	:to (new-gname 'YYPOS-1011-3) 
	:rate (s* (new-gname 'YYPOS-1010-3) (new-gname 'strainacc-mdr-yy)))


(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1001-1) 
	:to (new-gname 'YYPOS-1101-1) 
	:rate (s* (new-gname 'YYPOS-1001-1) (new-gname 'strainacc-inh-yy)))
(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1001-1) 
	:to (new-gname 'YYPOS-1011-1) 
	:rate (s* (new-gname 'YYPOS-1001-1) (new-gname 'strainacc-rif-yy)))
(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1001-1) 
	:to (new-gname 'YYPOS-1001-4) 
	:rate (s* (new-gname 'YYPOS-1001-1) (new-gname 'overgrowth-yy)))


(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-0110-2) 
	:to (new-gname 'YYPOS-0111-2) 
	:rate (s* (new-gname 'YYPOS-0110-2) (new-gname 'strainacc-mdr-yy)))


(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-0110-3) 
	:to (new-gname 'YYPOS-0111-3) 
	:rate (s* (new-gname 'YYPOS-0110-3) (new-gname 'strainacc-mdr-yy)))


(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-0101-2) 
	:to (new-gname 'YYPOS-0101-4) 
	:rate (s* (new-gname 'YYPOS-0101-2) (new-gname 'overgrowth-yy)))


(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-0011-3) 
	:to (new-gname 'YYPOS-0011-4) 
	:rate (s* (new-gname 'YYPOS-0011-3) (new-gname 'overgrowth-yy)))


(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1110-1) 
	:to (new-gname 'YYPOS-1110-2) 
	:rate (s* (new-gname 'YYPOS-1110-1) (new-gname 'overgrowth-yy) ))
(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1110-1) 
	:to (new-gname 'YYPOS-1110-3) 
	:rate (s* (new-gname 'YYPOS-1110-1) (new-gname 'overgrowth-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1110-2) 
	:to (new-gname 'YYPOS-1111-2) 
	:rate (s* (new-gname 'YYPOS-1110-2) (new-gname 'strainacc-mdr-yy)))


(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1110-3) 
	:to (new-gname 'YYPOS-1111-3) 
	:rate (s* (new-gname 'YYPOS-1110-3) (new-gname 'strainacc-mdr-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1101-1) 
	:to (new-gname 'YYPOS-1111-1) 
	:rate (s* (new-gname 'YYPOS-1101-1) (new-gname 'strainacc-rif-yy) ))
(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1101-1) 
	:to (new-gname 'YYPOS-1101-2) 
	:rate (s* (new-gname 'YYPOS-1101-1) (new-gname 'overgrowth-yy)))
(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1101-1) 
	:to (new-gname 'YYPOS-1101-4) 
	:rate (s* (new-gname 'YYPOS-1101-1) (new-gname 'overgrowth-yy)))


(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1101-2) 
	:to (new-gname 'YYPOS-1101-4) 
	:rate (s* (new-gname 'YYPOS-1101-2) (new-gname 'overgrowth-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1011-1) 
	:to (new-gname 'YYPOS-1111-1) 
	:rate (s* (new-gname 'YYPOS-1011-1) (new-gname 'strainacc-inh-yy)))
(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1011-1) 
	:to (new-gname 'YYPOS-1011-3) 
	:rate (s* (new-gname 'YYPOS-1011-1) (new-gname 'overgrowth-yy)))
(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1011-1) 
	:to (new-gname 'YYPOS-1011-4) 
	:rate (s* (new-gname 'YYPOS-1011-1) (new-gname 'overgrowth-yy)))


(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1011-3) 
	:to (new-gname 'YYPOS-1011-4) 
	:rate (s* (new-gname 'YYPOS-1011-3) (new-gname 'overgrowth-yy)))



(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-0111-2) 
	:to (new-gname 'YYPOS-0111-4) 
	:rate (s* (new-gname 'YYPOS-0111-2) (new-gname 'overgrowth-yy)))


(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-0111-3) 
	:to (new-gname 'YYPOS-0111-4) 
	:rate (s* (new-gname 'YYPOS-0111-3) (new-gname 'overgrowth-yy)))


(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1111-1) 
	:to (new-gname 'YYPOS-1111-2) 
	:rate (s* (new-gname 'YYPOS-1111-1) (new-gname 'overgrowth-yy)))
(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1111-1) 
	:to (new-gname 'YYPOS-1111-3) 
	:rate (s* (new-gname 'YYPOS-1111-1) (new-gname 'overgrowth-yy)))
(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1111-1) 
	:to (new-gname 'YYPOS-1111-4) 
	:rate (s* (new-gname 'YYPOS-1111-1) (new-gname 'overgrowth-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1111-2) 
	:to (new-gname 'YYPOS-1111-4) 
	:rate (s* (new-gname 'YYPOS-1111-2) (new-gname 'overgrowth-yy) ))

(add-transition tbmodel  
	 
	:from (new-gname 'YYPOS-1111-3) 
	:to (new-gname 'YYPOS-1111-4) 
	:rate (s* (new-gname 'YYPOS-1111-3) (new-gname 'overgrowth-yy) ))
@

YYNEG
<<>>=
(add-transition tbmodel  
	 
	:from  (new-gname 'YYNEG-1000-1) 
	:to (new-gname 'YYNEG-1100-1) 
	:rate (s* (new-gname 'YYNEG-1000-1) (new-gname 'strainacc-inh-yy)))
(add-transition tbmodel  
	 
	:from  (new-gname 'YYNEG-1000-1) 
	:to (new-gname 'YYNEG-1010-1) 
	:rate (s* (new-gname 'YYNEG-1000-1) (new-gname 'strainacc-rif-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-0100-2) 
	:to (new-gname 'YYNEG-0101-2) 
	:rate (s* (new-gname 'YYNEG-0100-2) (new-gname 'strainacc-mdr-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-0010-3) 
	:to (new-gname 'YYNEG-0011-3) 
	:rate (s* (new-gname 'YYNEG-0010-3) (new-gname 'strainacc-mdr-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1100-1) 
	:to (new-gname 'YYNEG-1110-1) 
	:rate (s* (new-gname 'YYNEG-1100-1) (new-gname 'strainacc-rif-yy)))
(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1100-1) 
	:to (new-gname 'YYNEG-1100-2) 
	:rate (s* (new-gname 'YYNEG-1100-1) (new-gname 'overgrowth-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1100-2) 
	:to (new-gname 'YYNEG-1101-2) 
	:rate (s* (new-gname 'YYNEG-1100-2) (new-gname 'strainacc-mdr-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1010-1) 
	:to (new-gname 'YYNEG-1110-1) 
	:rate (s* (new-gname 'YYNEG-1010-1) (new-gname 'strainacc-inh-yy)))
(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1010-1) 
	:to (new-gname 'YYNEG-1010-3) 
	:rate (s* (new-gname 'YYNEG-1010-1) (new-gname 'overgrowth-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1010-3) 
	:to (new-gname 'YYNEG-1011-3) 
	:rate (s* (new-gname 'YYNEG-1010-3) (new-gname 'strainacc-mdr-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1001-1) 
	:to (new-gname 'YYNEG-1101-1) 
	:rate (s* (new-gname 'YYNEG-1001-1) (new-gname 'strainacc-inh-yy)))
(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1001-1) 
	:to (new-gname 'YYNEG-1011-1) 
	:rate (s* (new-gname 'YYNEG-1001-1) (new-gname 'strainacc-rif-yy)))
(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1001-1) 
	:to (new-gname 'YYNEG-1001-4) 
	:rate (s* (new-gname 'YYNEG-1001-1) (new-gname 'overgrowth-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-0110-2) 
	:to (new-gname 'YYNEG-0111-2) 
	:rate (s* (new-gname 'YYNEG-0110-2) (new-gname 'strainacc-mdr-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-0110-3) 
	:to (new-gname 'YYNEG-0111-3) 
	:rate (s* (new-gname 'YYNEG-0110-3) (new-gname 'strainacc-mdr-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-0101-2) 
	:to (new-gname 'YYNEG-0101-4) 
	:rate (s* (new-gname 'YYNEG-0101-2) (new-gname 'overgrowth-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-0011-3) 
	:to (new-gname 'YYNEG-0011-4) 
	:rate (s* (new-gname 'YYNEG-0011-3) (new-gname 'overgrowth-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1110-1) 
	:to (new-gname 'YYNEG-1110-2) 
	:rate (s* (new-gname 'YYNEG-1110-1) (new-gname 'overgrowth-yy) ))
(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1110-1) 
	:to (new-gname 'YYNEG-1110-3) 
	:rate (s* (new-gname 'YYNEG-1110-1) (new-gname 'overgrowth-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1110-2) 
	:to (new-gname 'YYNEG-1111-2) 
	:rate (s* (new-gname 'YYNEG-1110-2) (new-gname 'strainacc-mdr-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1110-3) 
	:to (new-gname 'YYNEG-1111-3) 
	:rate (s* (new-gname 'YYNEG-1110-3) (new-gname 'strainacc-mdr-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1101-1) 
	:to (new-gname 'YYNEG-1111-1) 
	:rate (s* (new-gname 'YYNEG-1101-1) (new-gname 'strainacc-rif-yy) ))
(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1101-1) 
	:to (new-gname 'YYNEG-1101-2) 
	:rate (s* (new-gname 'YYNEG-1101-1) (new-gname 'overgrowth-yy)))
(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1101-1) 
	:to (new-gname 'YYNEG-1101-4) 
	:rate (s* (new-gname 'YYNEG-1101-1) (new-gname 'overgrowth-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1101-2) 
	:to (new-gname 'YYNEG-1101-4) 
	:rate (s* (new-gname 'YYNEG-1101-2) (new-gname 'overgrowth-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1011-1) 
	:to (new-gname 'YYNEG-1111-1) 
	:rate (s* (new-gname 'YYNEG-1011-1) (new-gname 'strainacc-inh-yy)))
(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1011-1) 
	:to (new-gname 'YYNEG-1011-3) 
	:rate (s* (new-gname 'YYNEG-1011-1) (new-gname 'overgrowth-yy)))
(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1011-1) 
	:to (new-gname 'YYNEG-1011-4) 
	:rate (s* (new-gname 'YYNEG-1011-1) (new-gname 'overgrowth-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1011-3) 
	:to (new-gname 'YYNEG-1011-4) 
	:rate (s* (new-gname 'YYNEG-1011-3) (new-gname 'overgrowth-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-0111-2) 
	:to (new-gname 'YYNEG-0111-4) 
	:rate (s* (new-gname 'YYNEG-0111-2) (new-gname 'overgrowth-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-0111-3) 
	:to (new-gname 'YYNEG-0111-4) 
	:rate (s* (new-gname 'YYNEG-0111-3) (new-gname 'overgrowth-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1111-1) 
	:to (new-gname 'YYNEG-1111-2) 
	:rate (s* (new-gname 'YYNEG-1111-1) (new-gname 'overgrowth-yy)))
(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1111-1) 
	:to (new-gname 'YYNEG-1111-3) 
	:rate (s* (new-gname 'YYNEG-1111-1) (new-gname 'overgrowth-yy)))
(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1111-1) 
	:to (new-gname 'YYNEG-1111-4) 
	:rate (s* (new-gname 'YYNEG-1111-1) (new-gname 'overgrowth-yy)))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1111-2) 
	:to (new-gname 'YYNEG-1111-4) 
	:rate (s* (new-gname 'YYNEG-1111-2) (new-gname 'overgrowth-yy) ))

(add-transition tbmodel  
	 
	:from (new-gname 'YYNEG-1111-3) 
	:to (new-gname 'YYNEG-1111-4) 
	:rate (s* (new-gname 'YYNEG-1111-3) (new-gname 'overgrowth-yy) ))
@

The next section adds in reinfection-related transitions 
<<>>=
(add-parameter tbmodel (new-gname 'lam1))
(add-parameter tbmodel (new-gname 'lam2))
(add-parameter tbmodel (new-gname 'lam3))
(add-parameter tbmodel (new-gname 'lam4))
(add-parameter tbmodel (new-gname 'pfast))
(add-parameter tbmodel (new-gname 'zeta))
@
Rationale:
LL-abcd-e denotes the number of people with a given infection;
  {\tt a} indicates if the person is infected with WT;
  {\tt b} indicates if the person is infected with INH-monoresistant;
  {\tt c} indicates if the person is infected with rifampin monoresistant;
  {\tt d} indicates if the person is infected with MDR TB;
  {\tt e} tells which strain is dominant: 1 WT, 2 INH mono, 3 RIF mono, 4 MDR.
The assumption is that if a person is in LL or LLPRIME, whatever strain is dominant is already causing
progression to disease.
Reinfection can serve (by our assumption) only to add new strains within the host
we assume that the reinfecting or superinfecting strain does not ever take over from the currently
predominating strain in a person in state LL or LLPRIME.
Reinfection by the same strain a person is already predominated by is a no-op 
We do assume individuals are partially protected from reinfection, but we do not separate this
probability by state, only previous infection status.
RR states behave like LL or LLPRIME.
 
For MM and MMPRIME, reinfection is modeled differently.  These individuals are NOT on the fast track to 
tuberculosis.  We make some simplifying assumptions as follows: 1. Reinfection could add a new strain
without changing which strain predominates and without moving the person into rapid progression, or 2. the
reinfecting strain could move the person into rapid progression and we assume that it becomes predominant.
It is understood that a reinfecting strain could in principle cause a quiescent infection to flare up, so
that the existing strain or previous dominant strain stayed dominant, but we choose to ignore this for now.
Here pfast is the chance the reinfection triggers fast progression, and we have the same zeta as before.

<<>>=
(add-transition tbmodel   :from (new-gname 'MM-1000-1) :to (new-gname 'LL-1000-1) :rate (s* (new-gname 'MM-1000-1) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1000-1) :to (new-gname 'MM-1100-1) :rate (s* (new-gname 'MM-1000-1) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1000-1) :to (new-gname 'LL-1100-2) :rate (s* (new-gname 'MM-1000-1) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1000-1) :to (new-gname 'MM-1010-1) :rate (s* (new-gname 'MM-1000-1) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1000-1) :to (new-gname 'LL-1010-3) :rate (s* (new-gname 'MM-1000-1) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1000-1) :to (new-gname 'MM-1001-1) :rate (s* (new-gname 'MM-1000-1) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1000-1) :to (new-gname 'LL-1001-4) :rate (s* (new-gname 'MM-1000-1) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-1000-1) :to (new-gname 'LL-1100-1) :rate (s* (new-gname 'LL-1000-1) (new-gname 'lam2) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LL-1000-1) :to (new-gname 'LL-1010-1) :rate (s* (new-gname 'LL-1000-1) (new-gname 'lam3) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LL-1000-1) :to (new-gname 'LL-1001-1) :rate (s* (new-gname 'LL-1000-1) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-0100-2) :to (new-gname 'MM-1100-2) :rate (s* (new-gname 'MM-0100-2) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0100-2) :to (new-gname 'LL-1100-1) :rate (s* (new-gname 'MM-0100-2) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0100-2) :to (new-gname 'LL-0100-2) :rate (s* (new-gname 'MM-0100-2) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0100-2) :to (new-gname 'MM-0110-2) :rate (s* (new-gname 'MM-0100-2) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0100-2) :to (new-gname 'LL-0110-3) :rate (s* (new-gname 'MM-0100-2) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0100-2) :to (new-gname 'MM-0101-2) :rate (s* (new-gname 'MM-0100-2) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0100-2) :to (new-gname 'LL-0101-4) :rate (s* (new-gname 'MM-0100-2) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-0100-2) :to (new-gname 'LL-1100-2) :rate (s* (new-gname 'LL-0100-2) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LL-0100-2) :to (new-gname 'LL-0110-2) :rate (s* (new-gname 'LL-0100-2) (new-gname 'lam3) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LL-0100-2) :to (new-gname 'LL-0101-2) :rate (s* (new-gname 'LL-0100-2) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-0010-3) :to (new-gname 'MM-1010-3) :rate (s* (new-gname 'MM-0010-3) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0010-3) :to (new-gname 'LL-1010-1) :rate (s* (new-gname 'MM-0010-3) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0010-3) :to (new-gname 'MM-0110-3) :rate (s* (new-gname 'MM-0010-3) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0010-3) :to (new-gname 'LL-0110-2) :rate (s* (new-gname 'MM-0010-3) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0010-3) :to (new-gname 'LL-0010-3) :rate (s* (new-gname 'MM-0010-3) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0010-3) :to (new-gname 'MM-0011-3) :rate (s* (new-gname 'MM-0010-3) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0010-3) :to (new-gname 'LL-0011-4) :rate (s* (new-gname 'MM-0010-3) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-0010-3) :to (new-gname 'LL-1010-3) :rate (s* (new-gname 'LL-0010-3) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LL-0010-3) :to (new-gname 'LL-0110-3) :rate (s* (new-gname 'LL-0010-3) (new-gname 'lam2) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LL-0010-3) :to (new-gname 'LL-0011-3) :rate (s* (new-gname 'LL-0010-3) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-0001-4) :to (new-gname 'MM-1001-4) :rate (s* (new-gname 'MM-0001-4) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0001-4) :to (new-gname 'LL-1001-1) :rate (s* (new-gname 'MM-0001-4) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0001-4) :to (new-gname 'MM-0101-4) :rate (s* (new-gname 'MM-0001-4) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0001-4) :to (new-gname 'LL-0101-2) :rate (s* (new-gname 'MM-0001-4) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0001-4) :to (new-gname 'MM-0011-4) :rate (s* (new-gname 'MM-0001-4) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0001-4) :to (new-gname 'LL-0011-3) :rate (s* (new-gname 'MM-0001-4) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0001-4) :to (new-gname 'LL-0001-4) :rate (s* (new-gname 'MM-0001-4) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-0001-4) :to (new-gname 'LL-1001-4) :rate (s* (new-gname 'LL-0001-4) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LL-0001-4) :to (new-gname 'LL-0101-4) :rate (s* (new-gname 'LL-0001-4) (new-gname 'lam2) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LL-0001-4) :to (new-gname 'LL-0011-4) :rate (s* (new-gname 'LL-0001-4) (new-gname 'lam3) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-1100-1) :to (new-gname 'LL-1100-1) :rate (s* (new-gname 'MM-1100-1) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1100-1) :to (new-gname 'LL-1100-2) :rate (s* (new-gname 'MM-1100-1) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1100-1) :to (new-gname 'MM-1110-1) :rate (s* (new-gname 'MM-1100-1) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1100-1) :to (new-gname 'LL-1110-3) :rate (s* (new-gname 'MM-1100-1) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1100-1) :to (new-gname 'MM-1101-1) :rate (s* (new-gname 'MM-1100-1) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1100-1) :to (new-gname 'LL-1101-4) :rate (s* (new-gname 'MM-1100-1) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-1100-1) :to (new-gname 'LL-1110-1) :rate (s* (new-gname 'LL-1100-1) (new-gname 'lam3) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LL-1100-1) :to (new-gname 'LL-1101-1) :rate (s* (new-gname 'LL-1100-1) (new-gname 'lam4) (new-gname 'zeta) ))

(add-transition tbmodel   :from (new-gname 'MM-1100-2) :to (new-gname 'LL-1100-1) :rate (s* (new-gname 'MM-1100-2) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1100-2) :to (new-gname 'LL-1100-2) :rate (s* (new-gname 'MM-1100-2) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1100-2) :to (new-gname 'MM-1110-2) :rate (s* (new-gname 'MM-1100-2) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1100-2) :to (new-gname 'LL-1110-3) :rate (s* (new-gname 'MM-1100-2) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1100-2) :to (new-gname 'MM-1101-2) :rate (s* (new-gname 'MM-1100-2) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1100-2) :to (new-gname 'LL-1101-4) :rate (s* (new-gname 'MM-1100-2) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-1100-2) :to (new-gname 'LL-1110-2) :rate (s* (new-gname 'LL-1100-2) (new-gname 'lam3) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LL-1100-2) :to (new-gname 'LL-1101-2) :rate (s* (new-gname 'LL-1100-2) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-1010-1) :to (new-gname 'LL-1010-1) :rate (s* (new-gname 'MM-1010-1) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1010-1) :to (new-gname 'MM-1110-1) :rate (s* (new-gname 'MM-1010-1) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1010-1) :to (new-gname 'LL-1110-2) :rate (s* (new-gname 'MM-1010-1) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1010-1) :to (new-gname 'LL-1010-3) :rate (s* (new-gname 'MM-1010-1) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1010-1) :to (new-gname 'MM-1011-1) :rate (s* (new-gname 'MM-1010-1) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1010-1) :to (new-gname 'LL-1011-4) :rate (s* (new-gname 'MM-1010-1) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-1010-1) :to (new-gname 'LL-1110-1) :rate (s* (new-gname 'LL-1010-1) (new-gname 'lam2) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LL-1010-1) :to (new-gname 'LL-1011-1) :rate (s* (new-gname 'LL-1010-1) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-1010-3) :to (new-gname 'LL-1010-1) :rate (s* (new-gname 'MM-1010-3) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1010-3) :to (new-gname 'MM-1110-3) :rate (s* (new-gname 'MM-1010-3) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1010-3) :to (new-gname 'LL-1110-2) :rate (s* (new-gname 'MM-1010-3) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1010-3) :to (new-gname 'LL-1010-3) :rate (s* (new-gname 'MM-1010-3) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1010-3) :to (new-gname 'MM-1011-3) :rate (s* (new-gname 'MM-1010-3) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1010-3) :to (new-gname 'LL-1011-4) :rate (s* (new-gname 'MM-1010-3) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-1010-3) :to (new-gname 'LL-1110-3) :rate (s* (new-gname 'LL-1010-3) (new-gname 'lam2) (new-gname 'zeta) ))
(add-transition tbmodel   :from (new-gname 'LL-1010-3) :to (new-gname 'LL-1011-3) :rate (s* (new-gname 'LL-1010-3) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-1001-1) :to (new-gname 'LL-1001-1) :rate (s* (new-gname 'MM-1001-1) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1001-1) :to (new-gname 'MM-1101-1) :rate (s* (new-gname 'MM-1001-1) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1001-1) :to (new-gname 'LL-1101-2) :rate (s* (new-gname 'MM-1001-1) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1001-1) :to (new-gname 'MM-1011-1) :rate (s* (new-gname 'MM-1001-1) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1001-1) :to (new-gname 'LL-1011-3) :rate (s* (new-gname 'MM-1001-1) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1001-1) :to (new-gname 'LL-1001-4) :rate (s* (new-gname 'MM-1001-1) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-1001-1) :to (new-gname 'LL-1101-1) :rate (s* (new-gname 'LL-1001-1) (new-gname 'lam2) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LL-1001-1) :to (new-gname 'LL-1011-1) :rate (s* (new-gname 'LL-1001-1) (new-gname 'lam3) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-1001-4) :to (new-gname 'LL-1001-1) :rate (s* (new-gname 'MM-1001-4) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1001-4) :to (new-gname 'MM-1101-4) :rate (s* (new-gname 'MM-1001-4) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1001-4) :to (new-gname 'LL-1101-2) :rate (s* (new-gname 'MM-1001-4) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1001-4) :to (new-gname 'MM-1011-4) :rate (s* (new-gname 'MM-1001-4) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1001-4) :to (new-gname 'LL-1011-3) :rate (s* (new-gname 'MM-1001-4) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1001-4) :to (new-gname 'LL-1001-4) :rate (s* (new-gname 'MM-1001-4) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-1001-4) :to (new-gname 'LL-1101-4) :rate (s* (new-gname 'LL-1001-4) (new-gname 'lam2) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LL-1001-4) :to (new-gname 'LL-1011-4) :rate (s* (new-gname 'LL-1001-4) (new-gname 'lam3) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-0110-2) :to (new-gname 'MM-1110-2) :rate (s* (new-gname 'MM-0110-2) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0110-2) :to (new-gname 'LL-1110-1) :rate (s* (new-gname 'MM-0110-2) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0110-2) :to (new-gname 'LL-0110-2) :rate (s* (new-gname 'MM-0110-2) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0110-2) :to (new-gname 'LL-0110-3) :rate (s* (new-gname 'MM-0110-2) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0110-2) :to (new-gname 'MM-0111-2) :rate (s* (new-gname 'MM-0110-2) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0110-2) :to (new-gname 'LL-0111-4) :rate (s* (new-gname 'MM-0110-2) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-0110-2) :to (new-gname 'LL-1110-2) :rate (s* (new-gname 'LL-0110-2) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LL-0110-2) :to (new-gname 'LL-0111-2) :rate (s* (new-gname 'LL-0110-2) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-0110-3) :to (new-gname 'MM-1110-3) :rate (s* (new-gname 'MM-0110-3) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0110-3) :to (new-gname 'LL-1110-1) :rate (s* (new-gname 'MM-0110-3) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0110-3) :to (new-gname 'LL-0110-2) :rate (s* (new-gname 'MM-0110-3) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0110-3) :to (new-gname 'LL-0110-3) :rate (s* (new-gname 'MM-0110-3) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0110-3) :to (new-gname 'MM-0111-3) :rate (s* (new-gname 'MM-0110-3) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0110-3) :to (new-gname 'LL-0111-4) :rate (s* (new-gname 'MM-0110-3) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-0110-3) :to (new-gname 'LL-1110-3) :rate (s* (new-gname 'LL-0110-3) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LL-0110-3) :to (new-gname 'LL-0111-3) :rate (s* (new-gname 'LL-0110-3) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-0101-2) :to (new-gname 'MM-1101-2) :rate (s* (new-gname 'MM-0101-2) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0101-2) :to (new-gname 'LL-1101-1) :rate (s* (new-gname 'MM-0101-2) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0101-2) :to (new-gname 'LL-0101-2) :rate (s* (new-gname 'MM-0101-2) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0101-2) :to (new-gname 'MM-0111-2) :rate (s* (new-gname 'MM-0101-2) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0101-2) :to (new-gname 'LL-0111-3) :rate (s* (new-gname 'MM-0101-2) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0101-2) :to (new-gname 'LL-0101-4) :rate (s* (new-gname 'MM-0101-2) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-0101-2) :to (new-gname 'LL-1101-2) :rate (s* (new-gname 'LL-0101-2) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LL-0101-2) :to (new-gname 'LL-0111-2) :rate (s* (new-gname 'LL-0101-2) (new-gname 'lam3) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-0101-4) :to (new-gname 'MM-1101-4) :rate (s* (new-gname 'MM-0101-4) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0101-4) :to (new-gname 'LL-1101-1) :rate (s* (new-gname 'MM-0101-4) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0101-4) :to (new-gname 'LL-0101-2) :rate (s* (new-gname 'MM-0101-4) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0101-4) :to (new-gname 'MM-0111-4) :rate (s* (new-gname 'MM-0101-4) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0101-4) :to (new-gname 'LL-0111-3) :rate (s* (new-gname 'MM-0101-4) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0101-4) :to (new-gname 'LL-0101-4) :rate (s* (new-gname 'MM-0101-4) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-0101-4) :to (new-gname 'LL-1101-4) :rate (s* (new-gname 'LL-0101-4) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LL-0101-4) :to (new-gname 'LL-0111-4) :rate (s* (new-gname 'LL-0101-4) (new-gname 'lam3) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-0011-3) :to (new-gname 'MM-1011-3) :rate (s* (new-gname 'MM-0011-3) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0011-3) :to (new-gname 'LL-1011-1) :rate (s* (new-gname 'MM-0011-3) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0011-3) :to (new-gname 'MM-0111-3) :rate (s* (new-gname 'MM-0011-3) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0011-3) :to (new-gname 'LL-0111-2) :rate (s* (new-gname 'MM-0011-3) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0011-3) :to (new-gname 'LL-0011-3) :rate (s* (new-gname 'MM-0011-3) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0011-3) :to (new-gname 'LL-0011-4) :rate (s* (new-gname 'MM-0011-3) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-0011-3) :to (new-gname 'LL-1011-3) :rate (s* (new-gname 'LL-0011-3) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LL-0011-3) :to (new-gname 'LL-0111-3) :rate (s* (new-gname 'LL-0011-3) (new-gname 'lam2) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-0011-4) :to (new-gname 'MM-1011-4) :rate (s* (new-gname 'MM-0011-4) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0011-4) :to (new-gname 'LL-1011-1) :rate (s* (new-gname 'MM-0011-4) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0011-4) :to (new-gname 'MM-0111-4) :rate (s* (new-gname 'MM-0011-4) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0011-4) :to (new-gname 'LL-0111-2) :rate (s* (new-gname 'MM-0011-4) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0011-4) :to (new-gname 'LL-0011-3) :rate (s* (new-gname 'MM-0011-4) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0011-4) :to (new-gname 'LL-0011-4) :rate (s* (new-gname 'MM-0011-4) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-0011-4) :to (new-gname 'LL-1011-4) :rate (s* (new-gname 'LL-0011-4) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LL-0011-4) :to (new-gname 'LL-0111-4) :rate (s* (new-gname 'LL-0011-4) (new-gname 'lam2) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-1110-1) :to (new-gname 'LL-1110-1) :rate (s* (new-gname 'MM-1110-1) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1110-1) :to (new-gname 'LL-1110-2) :rate (s* (new-gname 'MM-1110-1) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1110-1) :to (new-gname 'LL-1110-3) :rate (s* (new-gname 'MM-1110-1) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1110-1) :to (new-gname 'MM-1111-1) :rate (s* (new-gname 'MM-1110-1) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1110-1) :to (new-gname 'LL-1111-4) :rate (s* (new-gname 'MM-1110-1) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-1110-1) :to (new-gname 'LL-1111-1) :rate (s* (new-gname 'LL-1110-1) (new-gname 'lam4) (new-gname 'zeta) ))

(add-transition tbmodel   :from (new-gname 'MM-1110-2) :to (new-gname 'LL-1110-1) :rate (s* (new-gname 'MM-1110-2) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1110-2) :to (new-gname 'LL-1110-2) :rate (s* (new-gname 'MM-1110-2) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1110-2) :to (new-gname 'LL-1110-3) :rate (s* (new-gname 'MM-1110-2) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1110-2) :to (new-gname 'MM-1111-2) :rate (s* (new-gname 'MM-1110-2) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1110-2) :to (new-gname 'LL-1111-4) :rate (s* (new-gname 'MM-1110-2) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-1110-2) :to (new-gname 'LL-1111-2) :rate (s* (new-gname 'LL-1110-2) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-1110-3) :to (new-gname 'LL-1110-1) :rate (s* (new-gname 'MM-1110-3) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1110-3) :to (new-gname 'LL-1110-2) :rate (s* (new-gname 'MM-1110-3) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1110-3) :to (new-gname 'LL-1110-3) :rate (s* (new-gname 'MM-1110-3) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1110-3) :to (new-gname 'MM-1111-3) :rate (s* (new-gname 'MM-1110-3) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1110-3) :to (new-gname 'LL-1111-4) :rate (s* (new-gname 'MM-1110-3) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-1110-3) :to (new-gname 'LL-1111-3) :rate (s* (new-gname 'LL-1110-3) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-1101-1) :to (new-gname 'LL-1101-1) :rate (s* (new-gname 'MM-1101-1) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1101-1) :to (new-gname 'LL-1101-2) :rate (s* (new-gname 'MM-1101-1) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1101-1) :to (new-gname 'MM-1111-1) :rate (s* (new-gname 'MM-1101-1) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1101-1) :to (new-gname 'LL-1111-3) :rate (s* (new-gname 'MM-1101-1) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1101-1) :to (new-gname 'LL-1101-4) :rate (s* (new-gname 'MM-1101-1) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-1101-1) :to (new-gname 'LL-1111-1) :rate (s* (new-gname 'LL-1101-1) (new-gname 'lam3) (new-gname 'zeta) ))

(add-transition tbmodel   :from (new-gname 'MM-1101-2) :to (new-gname 'LL-1101-1) :rate (s* (new-gname 'MM-1101-2) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1101-2) :to (new-gname 'LL-1101-2) :rate (s* (new-gname 'MM-1101-2) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1101-2) :to (new-gname 'MM-1111-2) :rate (s* (new-gname 'MM-1101-2) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1101-2) :to (new-gname 'LL-1111-3) :rate (s* (new-gname 'MM-1101-2) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1101-2) :to (new-gname 'LL-1101-4) :rate (s* (new-gname 'MM-1101-2) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-1101-2) :to (new-gname 'LL-1111-2) :rate (s* (new-gname 'LL-1101-2) (new-gname 'lam3) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-1101-4) :to (new-gname 'LL-1101-1) :rate (s* (new-gname 'MM-1101-4) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1101-4) :to (new-gname 'LL-1101-2) :rate (s* (new-gname 'MM-1101-4) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1101-4) :to (new-gname 'MM-1111-4) :rate (s* (new-gname 'MM-1101-4) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1101-4) :to (new-gname 'LL-1111-3) :rate (s* (new-gname 'MM-1101-4) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1101-4) :to (new-gname 'LL-1101-4) :rate (s* (new-gname 'MM-1101-4) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-1101-4) :to (new-gname 'LL-1111-4) :rate (s* (new-gname 'LL-1101-4) (new-gname 'lam3) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-1011-1) :to (new-gname 'LL-1011-1) :rate (s* (new-gname 'MM-1011-1) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1011-1) :to (new-gname 'MM-1111-1) :rate (s* (new-gname 'MM-1011-1) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1011-1) :to (new-gname 'LL-1111-2) :rate (s* (new-gname 'MM-1011-1) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1011-1) :to (new-gname 'LL-1011-3) :rate (s* (new-gname 'MM-1011-1) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1011-1) :to (new-gname 'LL-1011-4) :rate (s* (new-gname 'MM-1011-1) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-1011-1) :to (new-gname 'LL-1111-1) :rate (s* (new-gname 'LL-1011-1) (new-gname 'lam2) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-1011-3) :to (new-gname 'LL-1011-1) :rate (s* (new-gname 'MM-1011-3) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1011-3) :to (new-gname 'MM-1111-3) :rate (s* (new-gname 'MM-1011-3) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1011-3) :to (new-gname 'LL-1111-2) :rate (s* (new-gname 'MM-1011-3) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1011-3) :to (new-gname 'LL-1011-3) :rate (s* (new-gname 'MM-1011-3) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1011-3) :to (new-gname 'LL-1011-4) :rate (s* (new-gname 'MM-1011-3) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-1011-3) :to (new-gname 'LL-1111-3) :rate (s* (new-gname 'LL-1011-3) (new-gname 'lam2) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-1011-4) :to (new-gname 'LL-1011-1) :rate (s* (new-gname 'MM-1011-4) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1011-4) :to (new-gname 'MM-1111-4) :rate (s* (new-gname 'MM-1011-4) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-1011-4) :to (new-gname 'LL-1111-2) :rate (s* (new-gname 'MM-1011-4) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1011-4) :to (new-gname 'LL-1011-3) :rate (s* (new-gname 'MM-1011-4) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1011-4) :to (new-gname 'LL-1011-4) :rate (s* (new-gname 'MM-1011-4) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-1011-4) :to (new-gname 'LL-1111-4) :rate (s* (new-gname 'LL-1011-4) (new-gname 'lam2) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-0111-2) :to (new-gname 'MM-1111-2) :rate (s* (new-gname 'MM-0111-2) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0111-2) :to (new-gname 'LL-1111-1) :rate (s* (new-gname 'MM-0111-2) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0111-2) :to (new-gname 'LL-0111-2) :rate (s* (new-gname 'MM-0111-2) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0111-2) :to (new-gname 'LL-0111-3) :rate (s* (new-gname 'MM-0111-2) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0111-2) :to (new-gname 'LL-0111-4) :rate (s* (new-gname 'MM-0111-2) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-0111-2) :to (new-gname 'LL-1111-2) :rate (s* (new-gname 'LL-0111-2) (new-gname 'lam1) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-0111-3) :to (new-gname 'MM-1111-3) :rate (s* (new-gname 'MM-0111-3) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0111-3) :to (new-gname 'LL-1111-1) :rate (s* (new-gname 'MM-0111-3) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0111-3) :to (new-gname 'LL-0111-2) :rate (s* (new-gname 'MM-0111-3) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0111-3) :to (new-gname 'LL-0111-3) :rate (s* (new-gname 'MM-0111-3) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0111-3) :to (new-gname 'LL-0111-4) :rate (s* (new-gname 'MM-0111-3) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-0111-3) :to (new-gname 'LL-1111-3) :rate (s* (new-gname 'LL-0111-3) (new-gname 'lam1) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-0111-4) :to (new-gname 'MM-1111-4) :rate (s* (new-gname 'MM-0111-4) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MM-0111-4) :to (new-gname 'LL-1111-1) :rate (s* (new-gname 'MM-0111-4) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0111-4) :to (new-gname 'LL-0111-2) :rate (s* (new-gname 'MM-0111-4) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0111-4) :to (new-gname 'LL-0111-3) :rate (s* (new-gname 'MM-0111-4) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-0111-4) :to (new-gname 'LL-0111-4) :rate (s* (new-gname 'MM-0111-4) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LL-0111-4) :to (new-gname 'LL-1111-4) :rate (s* (new-gname 'LL-0111-4) (new-gname 'lam1) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MM-1111-1) :to (new-gname 'LL-1111-1) :rate (s* (new-gname 'MM-1111-1) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1111-1) :to (new-gname 'LL-1111-2) :rate (s* (new-gname 'MM-1111-1) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1111-1) :to (new-gname 'LL-1111-3) :rate (s* (new-gname 'MM-1111-1) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1111-1) :to (new-gname 'LL-1111-4) :rate (s* (new-gname 'MM-1111-1) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'MM-1111-2) :to (new-gname 'LL-1111-1) :rate (s* (new-gname 'MM-1111-2) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1111-2) :to (new-gname 'LL-1111-2) :rate (s* (new-gname 'MM-1111-2) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1111-2) :to (new-gname 'LL-1111-3) :rate (s* (new-gname 'MM-1111-2) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1111-2) :to (new-gname 'LL-1111-4) :rate (s* (new-gname 'MM-1111-2) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'MM-1111-3) :to (new-gname 'LL-1111-1) :rate (s* (new-gname 'MM-1111-3) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1111-3) :to (new-gname 'LL-1111-2) :rate (s* (new-gname 'MM-1111-3) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1111-3) :to (new-gname 'LL-1111-3) :rate (s* (new-gname 'MM-1111-3) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1111-3) :to (new-gname 'LL-1111-4) :rate (s* (new-gname 'MM-1111-3) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'MM-1111-4) :to (new-gname 'LL-1111-1) :rate (s* (new-gname 'MM-1111-4) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1111-4) :to (new-gname 'LL-1111-2) :rate (s* (new-gname 'MM-1111-4) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1111-4) :to (new-gname 'LL-1111-3) :rate (s* (new-gname 'MM-1111-4) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MM-1111-4) :to (new-gname 'LL-1111-4) :rate (s* (new-gname 'MM-1111-4) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))
@

The exact same transitions are replicated for the LLPRIME and MMPRIME states
<<>>=
(add-transition tbmodel   :from (new-gname 'MMPRIME-1000-1) :to (new-gname 'LLPRIME-1000-1) :rate (s* (new-gname 'MMPRIME-1000-1) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1000-1) :to (new-gname 'MMPRIME-1100-1) :rate (s* (new-gname 'MMPRIME-1000-1) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1000-1) :to (new-gname 'LLPRIME-1100-2) :rate (s* (new-gname 'MMPRIME-1000-1) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1000-1) :to (new-gname 'MMPRIME-1010-1) :rate (s* (new-gname 'MMPRIME-1000-1) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1000-1) :to (new-gname 'LLPRIME-1010-3) :rate (s* (new-gname 'MMPRIME-1000-1) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1000-1) :to (new-gname 'MMPRIME-1001-1) :rate (s* (new-gname 'MMPRIME-1000-1) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1000-1) :to (new-gname 'LLPRIME-1001-4) :rate (s* (new-gname 'MMPRIME-1000-1) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-1000-1) :to (new-gname 'LLPRIME-1100-1) :rate (s* (new-gname 'LLPRIME-1000-1) (new-gname 'lam2) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LLPRIME-1000-1) :to (new-gname 'LLPRIME-1010-1) :rate (s* (new-gname 'LLPRIME-1000-1) (new-gname 'lam3) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LLPRIME-1000-1) :to (new-gname 'LLPRIME-1001-1) :rate (s* (new-gname 'LLPRIME-1000-1) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-0100-2) :to (new-gname 'MMPRIME-1100-2) :rate (s* (new-gname 'MMPRIME-0100-2) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0100-2) :to (new-gname 'LLPRIME-1100-1) :rate (s* (new-gname 'MMPRIME-0100-2) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0100-2) :to (new-gname 'LLPRIME-0100-2) :rate (s* (new-gname 'MMPRIME-0100-2) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0100-2) :to (new-gname 'MMPRIME-0110-2) :rate (s* (new-gname 'MMPRIME-0100-2) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0100-2) :to (new-gname 'LLPRIME-0110-3) :rate (s* (new-gname 'MMPRIME-0100-2) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0100-2) :to (new-gname 'MMPRIME-0101-2) :rate (s* (new-gname 'MMPRIME-0100-2) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0100-2) :to (new-gname 'LLPRIME-0101-4) :rate (s* (new-gname 'MMPRIME-0100-2) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-0100-2) :to (new-gname 'LLPRIME-1100-2) :rate (s* (new-gname 'LLPRIME-0100-2) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LLPRIME-0100-2) :to (new-gname 'LLPRIME-0110-2) :rate (s* (new-gname 'LLPRIME-0100-2) (new-gname 'lam3) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LLPRIME-0100-2) :to (new-gname 'LLPRIME-0101-2) :rate (s* (new-gname 'LLPRIME-0100-2) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-0010-3) :to (new-gname 'MMPRIME-1010-3) :rate (s* (new-gname 'MMPRIME-0010-3) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0010-3) :to (new-gname 'LLPRIME-1010-1) :rate (s* (new-gname 'MMPRIME-0010-3) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0010-3) :to (new-gname 'MMPRIME-0110-3) :rate (s* (new-gname 'MMPRIME-0010-3) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0010-3) :to (new-gname 'LLPRIME-0110-2) :rate (s* (new-gname 'MMPRIME-0010-3) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0010-3) :to (new-gname 'LLPRIME-0010-3) :rate (s* (new-gname 'MMPRIME-0010-3) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0010-3) :to (new-gname 'MMPRIME-0011-3) :rate (s* (new-gname 'MMPRIME-0010-3) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0010-3) :to (new-gname 'LLPRIME-0011-4) :rate (s* (new-gname 'MMPRIME-0010-3) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-0010-3) :to (new-gname 'LLPRIME-1010-3) :rate (s* (new-gname 'LLPRIME-0010-3) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LLPRIME-0010-3) :to (new-gname 'LLPRIME-0110-3) :rate (s* (new-gname 'LLPRIME-0010-3) (new-gname 'lam2) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LLPRIME-0010-3) :to (new-gname 'LLPRIME-0011-3) :rate (s* (new-gname 'LLPRIME-0010-3) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-0001-4) :to (new-gname 'MMPRIME-1001-4) :rate (s* (new-gname 'MMPRIME-0001-4) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0001-4) :to (new-gname 'LLPRIME-1001-1) :rate (s* (new-gname 'MMPRIME-0001-4) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0001-4) :to (new-gname 'MMPRIME-0101-4) :rate (s* (new-gname 'MMPRIME-0001-4) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0001-4) :to (new-gname 'LLPRIME-0101-2) :rate (s* (new-gname 'MMPRIME-0001-4) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0001-4) :to (new-gname 'MMPRIME-0011-4) :rate (s* (new-gname 'MMPRIME-0001-4) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0001-4) :to (new-gname 'LLPRIME-0011-3) :rate (s* (new-gname 'MMPRIME-0001-4) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0001-4) :to (new-gname 'LLPRIME-0001-4) :rate (s* (new-gname 'MMPRIME-0001-4) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-0001-4) :to (new-gname 'LLPRIME-1001-4) :rate (s* (new-gname 'LLPRIME-0001-4) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LLPRIME-0001-4) :to (new-gname 'LLPRIME-0101-4) :rate (s* (new-gname 'LLPRIME-0001-4) (new-gname 'lam2) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LLPRIME-0001-4) :to (new-gname 'LLPRIME-0011-4) :rate (s* (new-gname 'LLPRIME-0001-4) (new-gname 'lam3) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-1100-1) :to (new-gname 'LLPRIME-1100-1) :rate (s* (new-gname 'MMPRIME-1100-1) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1100-1) :to (new-gname 'LLPRIME-1100-2) :rate (s* (new-gname 'MMPRIME-1100-1) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1100-1) :to (new-gname 'MMPRIME-1110-1) :rate (s* (new-gname 'MMPRIME-1100-1) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1100-1) :to (new-gname 'LLPRIME-1110-3) :rate (s* (new-gname 'MMPRIME-1100-1) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1100-1) :to (new-gname 'MMPRIME-1101-1) :rate (s* (new-gname 'MMPRIME-1100-1) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1100-1) :to (new-gname 'LLPRIME-1101-4) :rate (s* (new-gname 'MMPRIME-1100-1) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-1100-1) :to (new-gname 'LLPRIME-1110-1) :rate (s* (new-gname 'LLPRIME-1100-1) (new-gname 'lam3) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LLPRIME-1100-1) :to (new-gname 'LLPRIME-1101-1) :rate (s* (new-gname 'LLPRIME-1100-1) (new-gname 'lam4) (new-gname 'zeta) ))

(add-transition tbmodel   :from (new-gname 'MMPRIME-1100-2) :to (new-gname 'LLPRIME-1100-1) :rate (s* (new-gname 'MMPRIME-1100-2) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1100-2) :to (new-gname 'LLPRIME-1100-2) :rate (s* (new-gname 'MMPRIME-1100-2) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1100-2) :to (new-gname 'MMPRIME-1110-2) :rate (s* (new-gname 'MMPRIME-1100-2) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1100-2) :to (new-gname 'LLPRIME-1110-3) :rate (s* (new-gname 'MMPRIME-1100-2) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1100-2) :to (new-gname 'MMPRIME-1101-2) :rate (s* (new-gname 'MMPRIME-1100-2) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1100-2) :to (new-gname 'LLPRIME-1101-4) :rate (s* (new-gname 'MMPRIME-1100-2) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-1100-2) :to (new-gname 'LLPRIME-1110-2) :rate (s* (new-gname 'LLPRIME-1100-2) (new-gname 'lam3) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LLPRIME-1100-2) :to (new-gname 'LLPRIME-1101-2) :rate (s* (new-gname 'LLPRIME-1100-2) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-1010-1) :to (new-gname 'LLPRIME-1010-1) :rate (s* (new-gname 'MMPRIME-1010-1) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1010-1) :to (new-gname 'MMPRIME-1110-1) :rate (s* (new-gname 'MMPRIME-1010-1) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1010-1) :to (new-gname 'LLPRIME-1110-2) :rate (s* (new-gname 'MMPRIME-1010-1) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1010-1) :to (new-gname 'LLPRIME-1010-3) :rate (s* (new-gname 'MMPRIME-1010-1) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1010-1) :to (new-gname 'MMPRIME-1011-1) :rate (s* (new-gname 'MMPRIME-1010-1) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1010-1) :to (new-gname 'LLPRIME-1011-4) :rate (s* (new-gname 'MMPRIME-1010-1) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-1010-1) :to (new-gname 'LLPRIME-1110-1) :rate (s* (new-gname 'LLPRIME-1010-1) (new-gname 'lam2) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LLPRIME-1010-1) :to (new-gname 'LLPRIME-1011-1) :rate (s* (new-gname 'LLPRIME-1010-1) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-1010-3) :to (new-gname 'LLPRIME-1010-1) :rate (s* (new-gname 'MMPRIME-1010-3) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1010-3) :to (new-gname 'MMPRIME-1110-3) :rate (s* (new-gname 'MMPRIME-1010-3) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1010-3) :to (new-gname 'LLPRIME-1110-2) :rate (s* (new-gname 'MMPRIME-1010-3) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1010-3) :to (new-gname 'LLPRIME-1010-3) :rate (s* (new-gname 'MMPRIME-1010-3) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1010-3) :to (new-gname 'MMPRIME-1011-3) :rate (s* (new-gname 'MMPRIME-1010-3) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1010-3) :to (new-gname 'LLPRIME-1011-4) :rate (s* (new-gname 'MMPRIME-1010-3) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-1010-3) :to (new-gname 'LLPRIME-1110-3) :rate (s* (new-gname 'LLPRIME-1010-3) (new-gname 'lam2) (new-gname 'zeta) ))
(add-transition tbmodel   :from (new-gname 'LLPRIME-1010-3) :to (new-gname 'LLPRIME-1011-3) :rate (s* (new-gname 'LLPRIME-1010-3) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-1001-1) :to (new-gname 'LLPRIME-1001-1) :rate (s* (new-gname 'MMPRIME-1001-1) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1001-1) :to (new-gname 'MMPRIME-1101-1) :rate (s* (new-gname 'MMPRIME-1001-1) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1001-1) :to (new-gname 'LLPRIME-1101-2) :rate (s* (new-gname 'MMPRIME-1001-1) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1001-1) :to (new-gname 'MMPRIME-1011-1) :rate (s* (new-gname 'MMPRIME-1001-1) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1001-1) :to (new-gname 'LLPRIME-1011-3) :rate (s* (new-gname 'MMPRIME-1001-1) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1001-1) :to (new-gname 'LLPRIME-1001-4) :rate (s* (new-gname 'MMPRIME-1001-1) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-1001-1) :to (new-gname 'LLPRIME-1101-1) :rate (s* (new-gname 'LLPRIME-1001-1) (new-gname 'lam2) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LLPRIME-1001-1) :to (new-gname 'LLPRIME-1011-1) :rate (s* (new-gname 'LLPRIME-1001-1) (new-gname 'lam3) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-1001-4) :to (new-gname 'LLPRIME-1001-1) :rate (s* (new-gname 'MMPRIME-1001-4) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1001-4) :to (new-gname 'MMPRIME-1101-4) :rate (s* (new-gname 'MMPRIME-1001-4) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1001-4) :to (new-gname 'LLPRIME-1101-2) :rate (s* (new-gname 'MMPRIME-1001-4) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1001-4) :to (new-gname 'MMPRIME-1011-4) :rate (s* (new-gname 'MMPRIME-1001-4) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1001-4) :to (new-gname 'LLPRIME-1011-3) :rate (s* (new-gname 'MMPRIME-1001-4) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1001-4) :to (new-gname 'LLPRIME-1001-4) :rate (s* (new-gname 'MMPRIME-1001-4) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-1001-4) :to (new-gname 'LLPRIME-1101-4) :rate (s* (new-gname 'LLPRIME-1001-4) (new-gname 'lam2) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LLPRIME-1001-4) :to (new-gname 'LLPRIME-1011-4) :rate (s* (new-gname 'LLPRIME-1001-4) (new-gname 'lam3) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-0110-2) :to (new-gname 'MMPRIME-1110-2) :rate (s* (new-gname 'MMPRIME-0110-2) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0110-2) :to (new-gname 'LLPRIME-1110-1) :rate (s* (new-gname 'MMPRIME-0110-2) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0110-2) :to (new-gname 'LLPRIME-0110-2) :rate (s* (new-gname 'MMPRIME-0110-2) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0110-2) :to (new-gname 'LLPRIME-0110-3) :rate (s* (new-gname 'MMPRIME-0110-2) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0110-2) :to (new-gname 'MMPRIME-0111-2) :rate (s* (new-gname 'MMPRIME-0110-2) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0110-2) :to (new-gname 'LLPRIME-0111-4) :rate (s* (new-gname 'MMPRIME-0110-2) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-0110-2) :to (new-gname 'LLPRIME-1110-2) :rate (s* (new-gname 'LLPRIME-0110-2) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LLPRIME-0110-2) :to (new-gname 'LLPRIME-0111-2) :rate (s* (new-gname 'LLPRIME-0110-2) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-0110-3) :to (new-gname 'MMPRIME-1110-3) :rate (s* (new-gname 'MMPRIME-0110-3) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0110-3) :to (new-gname 'LLPRIME-1110-1) :rate (s* (new-gname 'MMPRIME-0110-3) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0110-3) :to (new-gname 'LLPRIME-0110-2) :rate (s* (new-gname 'MMPRIME-0110-3) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0110-3) :to (new-gname 'LLPRIME-0110-3) :rate (s* (new-gname 'MMPRIME-0110-3) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0110-3) :to (new-gname 'MMPRIME-0111-3) :rate (s* (new-gname 'MMPRIME-0110-3) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0110-3) :to (new-gname 'LLPRIME-0111-4) :rate (s* (new-gname 'MMPRIME-0110-3) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-0110-3) :to (new-gname 'LLPRIME-1110-3) :rate (s* (new-gname 'LLPRIME-0110-3) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LLPRIME-0110-3) :to (new-gname 'LLPRIME-0111-3) :rate (s* (new-gname 'LLPRIME-0110-3) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-0101-2) :to (new-gname 'MMPRIME-1101-2) :rate (s* (new-gname 'MMPRIME-0101-2) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0101-2) :to (new-gname 'LLPRIME-1101-1) :rate (s* (new-gname 'MMPRIME-0101-2) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0101-2) :to (new-gname 'LLPRIME-0101-2) :rate (s* (new-gname 'MMPRIME-0101-2) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0101-2) :to (new-gname 'MMPRIME-0111-2) :rate (s* (new-gname 'MMPRIME-0101-2) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0101-2) :to (new-gname 'LLPRIME-0111-3) :rate (s* (new-gname 'MMPRIME-0101-2) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0101-2) :to (new-gname 'LLPRIME-0101-4) :rate (s* (new-gname 'MMPRIME-0101-2) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-0101-2) :to (new-gname 'LLPRIME-1101-2) :rate (s* (new-gname 'LLPRIME-0101-2) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LLPRIME-0101-2) :to (new-gname 'LLPRIME-0111-2) :rate (s* (new-gname 'LLPRIME-0101-2) (new-gname 'lam3) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-0101-4) :to (new-gname 'MMPRIME-1101-4) :rate (s* (new-gname 'MMPRIME-0101-4) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0101-4) :to (new-gname 'LLPRIME-1101-1) :rate (s* (new-gname 'MMPRIME-0101-4) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0101-4) :to (new-gname 'LLPRIME-0101-2) :rate (s* (new-gname 'MMPRIME-0101-4) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0101-4) :to (new-gname 'MMPRIME-0111-4) :rate (s* (new-gname 'MMPRIME-0101-4) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0101-4) :to (new-gname 'LLPRIME-0111-3) :rate (s* (new-gname 'MMPRIME-0101-4) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0101-4) :to (new-gname 'LLPRIME-0101-4) :rate (s* (new-gname 'MMPRIME-0101-4) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-0101-4) :to (new-gname 'LLPRIME-1101-4) :rate (s* (new-gname 'LLPRIME-0101-4) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LLPRIME-0101-4) :to (new-gname 'LLPRIME-0111-4) :rate (s* (new-gname 'LLPRIME-0101-4) (new-gname 'lam3) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-0011-3) :to (new-gname 'MMPRIME-1011-3) :rate (s* (new-gname 'MMPRIME-0011-3) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0011-3) :to (new-gname 'LLPRIME-1011-1) :rate (s* (new-gname 'MMPRIME-0011-3) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0011-3) :to (new-gname 'MMPRIME-0111-3) :rate (s* (new-gname 'MMPRIME-0011-3) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0011-3) :to (new-gname 'LLPRIME-0111-2) :rate (s* (new-gname 'MMPRIME-0011-3) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0011-3) :to (new-gname 'LLPRIME-0011-3) :rate (s* (new-gname 'MMPRIME-0011-3) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0011-3) :to (new-gname 'LLPRIME-0011-4) :rate (s* (new-gname 'MMPRIME-0011-3) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-0011-3) :to (new-gname 'LLPRIME-1011-3) :rate (s* (new-gname 'LLPRIME-0011-3) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LLPRIME-0011-3) :to (new-gname 'LLPRIME-0111-3) :rate (s* (new-gname 'LLPRIME-0011-3) (new-gname 'lam2) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-0011-4) :to (new-gname 'MMPRIME-1011-4) :rate (s* (new-gname 'MMPRIME-0011-4) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0011-4) :to (new-gname 'LLPRIME-1011-1) :rate (s* (new-gname 'MMPRIME-0011-4) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0011-4) :to (new-gname 'MMPRIME-0111-4) :rate (s* (new-gname 'MMPRIME-0011-4) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0011-4) :to (new-gname 'LLPRIME-0111-2) :rate (s* (new-gname 'MMPRIME-0011-4) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0011-4) :to (new-gname 'LLPRIME-0011-3) :rate (s* (new-gname 'MMPRIME-0011-4) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0011-4) :to (new-gname 'LLPRIME-0011-4) :rate (s* (new-gname 'MMPRIME-0011-4) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-0011-4) :to (new-gname 'LLPRIME-1011-4) :rate (s* (new-gname 'LLPRIME-0011-4) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'LLPRIME-0011-4) :to (new-gname 'LLPRIME-0111-4) :rate (s* (new-gname 'LLPRIME-0011-4) (new-gname 'lam2) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-1110-1) :to (new-gname 'LLPRIME-1110-1) :rate (s* (new-gname 'MMPRIME-1110-1) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1110-1) :to (new-gname 'LLPRIME-1110-2) :rate (s* (new-gname 'MMPRIME-1110-1) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1110-1) :to (new-gname 'LLPRIME-1110-3) :rate (s* (new-gname 'MMPRIME-1110-1) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1110-1) :to (new-gname 'MMPRIME-1111-1) :rate (s* (new-gname 'MMPRIME-1110-1) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1110-1) :to (new-gname 'LLPRIME-1111-4) :rate (s* (new-gname 'MMPRIME-1110-1) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-1110-1) :to (new-gname 'LLPRIME-1111-1) :rate (s* (new-gname 'LLPRIME-1110-1) (new-gname 'lam4) (new-gname 'zeta) ))

(add-transition tbmodel   :from (new-gname 'MMPRIME-1110-2) :to (new-gname 'LLPRIME-1110-1) :rate (s* (new-gname 'MMPRIME-1110-2) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1110-2) :to (new-gname 'LLPRIME-1110-2) :rate (s* (new-gname 'MMPRIME-1110-2) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1110-2) :to (new-gname 'LLPRIME-1110-3) :rate (s* (new-gname 'MMPRIME-1110-2) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1110-2) :to (new-gname 'MMPRIME-1111-2) :rate (s* (new-gname 'MMPRIME-1110-2) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1110-2) :to (new-gname 'LLPRIME-1111-4) :rate (s* (new-gname 'MMPRIME-1110-2) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-1110-2) :to (new-gname 'LLPRIME-1111-2) :rate (s* (new-gname 'LLPRIME-1110-2) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-1110-3) :to (new-gname 'LLPRIME-1110-1) :rate (s* (new-gname 'MMPRIME-1110-3) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1110-3) :to (new-gname 'LLPRIME-1110-2) :rate (s* (new-gname 'MMPRIME-1110-3) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1110-3) :to (new-gname 'LLPRIME-1110-3) :rate (s* (new-gname 'MMPRIME-1110-3) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1110-3) :to (new-gname 'MMPRIME-1111-3) :rate (s* (new-gname 'MMPRIME-1110-3) (new-gname 'lam4) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1110-3) :to (new-gname 'LLPRIME-1111-4) :rate (s* (new-gname 'MMPRIME-1110-3) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-1110-3) :to (new-gname 'LLPRIME-1111-3) :rate (s* (new-gname 'LLPRIME-1110-3) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-1101-1) :to (new-gname 'LLPRIME-1101-1) :rate (s* (new-gname 'MMPRIME-1101-1) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1101-1) :to (new-gname 'LLPRIME-1101-2) :rate (s* (new-gname 'MMPRIME-1101-1) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1101-1) :to (new-gname 'MMPRIME-1111-1) :rate (s* (new-gname 'MMPRIME-1101-1) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1101-1) :to (new-gname 'LLPRIME-1111-3) :rate (s* (new-gname 'MMPRIME-1101-1) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1101-1) :to (new-gname 'LLPRIME-1101-4) :rate (s* (new-gname 'MMPRIME-1101-1) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-1101-1) :to (new-gname 'LLPRIME-1111-1) :rate (s* (new-gname 'LLPRIME-1101-1) (new-gname 'lam3) (new-gname 'zeta) ))

(add-transition tbmodel   :from (new-gname 'MMPRIME-1101-2) :to (new-gname 'LLPRIME-1101-1) :rate (s* (new-gname 'MMPRIME-1101-2) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1101-2) :to (new-gname 'LLPRIME-1101-2) :rate (s* (new-gname 'MMPRIME-1101-2) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1101-2) :to (new-gname 'MMPRIME-1111-2) :rate (s* (new-gname 'MMPRIME-1101-2) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1101-2) :to (new-gname 'LLPRIME-1111-3) :rate (s* (new-gname 'MMPRIME-1101-2) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1101-2) :to (new-gname 'LLPRIME-1101-4) :rate (s* (new-gname 'MMPRIME-1101-2) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-1101-2) :to (new-gname 'LLPRIME-1111-2) :rate (s* (new-gname 'LLPRIME-1101-2) (new-gname 'lam3) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-1101-4) :to (new-gname 'LLPRIME-1101-1) :rate (s* (new-gname 'MMPRIME-1101-4) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1101-4) :to (new-gname 'LLPRIME-1101-2) :rate (s* (new-gname 'MMPRIME-1101-4) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1101-4) :to (new-gname 'MMPRIME-1111-4) :rate (s* (new-gname 'MMPRIME-1101-4) (new-gname 'lam3) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1101-4) :to (new-gname 'LLPRIME-1111-3) :rate (s* (new-gname 'MMPRIME-1101-4) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1101-4) :to (new-gname 'LLPRIME-1101-4) :rate (s* (new-gname 'MMPRIME-1101-4) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-1101-4) :to (new-gname 'LLPRIME-1111-4) :rate (s* (new-gname 'LLPRIME-1101-4) (new-gname 'lam3) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-1011-1) :to (new-gname 'LLPRIME-1011-1) :rate (s* (new-gname 'MMPRIME-1011-1) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1011-1) :to (new-gname 'MMPRIME-1111-1) :rate (s* (new-gname 'MMPRIME-1011-1) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1011-1) :to (new-gname 'LLPRIME-1111-2) :rate (s* (new-gname 'MMPRIME-1011-1) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1011-1) :to (new-gname 'LLPRIME-1011-3) :rate (s* (new-gname 'MMPRIME-1011-1) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1011-1) :to (new-gname 'LLPRIME-1011-4) :rate (s* (new-gname 'MMPRIME-1011-1) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-1011-1) :to (new-gname 'LLPRIME-1111-1) :rate (s* (new-gname 'LLPRIME-1011-1) (new-gname 'lam2) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-1011-3) :to (new-gname 'LLPRIME-1011-1) :rate (s* (new-gname 'MMPRIME-1011-3) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1011-3) :to (new-gname 'MMPRIME-1111-3) :rate (s* (new-gname 'MMPRIME-1011-3) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1011-3) :to (new-gname 'LLPRIME-1111-2) :rate (s* (new-gname 'MMPRIME-1011-3) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1011-3) :to (new-gname 'LLPRIME-1011-3) :rate (s* (new-gname 'MMPRIME-1011-3) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1011-3) :to (new-gname 'LLPRIME-1011-4) :rate (s* (new-gname 'MMPRIME-1011-3) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-1011-3) :to (new-gname 'LLPRIME-1111-3) :rate (s* (new-gname 'LLPRIME-1011-3) (new-gname 'lam2) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-1011-4) :to (new-gname 'LLPRIME-1011-1) :rate (s* (new-gname 'MMPRIME-1011-4) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1011-4) :to (new-gname 'MMPRIME-1111-4) :rate (s* (new-gname 'MMPRIME-1011-4) (new-gname 'lam2) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1011-4) :to (new-gname 'LLPRIME-1111-2) :rate (s* (new-gname 'MMPRIME-1011-4) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1011-4) :to (new-gname 'LLPRIME-1011-3) :rate (s* (new-gname 'MMPRIME-1011-4) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1011-4) :to (new-gname 'LLPRIME-1011-4) :rate (s* (new-gname 'MMPRIME-1011-4) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-1011-4) :to (new-gname 'LLPRIME-1111-4) :rate (s* (new-gname 'LLPRIME-1011-4) (new-gname 'lam2) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-0111-2) :to (new-gname 'MMPRIME-1111-2) :rate (s* (new-gname 'MMPRIME-0111-2) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0111-2) :to (new-gname 'LLPRIME-1111-1) :rate (s* (new-gname 'MMPRIME-0111-2) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0111-2) :to (new-gname 'LLPRIME-0111-2) :rate (s* (new-gname 'MMPRIME-0111-2) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0111-2) :to (new-gname 'LLPRIME-0111-3) :rate (s* (new-gname 'MMPRIME-0111-2) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0111-2) :to (new-gname 'LLPRIME-0111-4) :rate (s* (new-gname 'MMPRIME-0111-2) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-0111-2) :to (new-gname 'LLPRIME-1111-2) :rate (s* (new-gname 'LLPRIME-0111-2) (new-gname 'lam1) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-0111-3) :to (new-gname 'MMPRIME-1111-3) :rate (s* (new-gname 'MMPRIME-0111-3) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0111-3) :to (new-gname 'LLPRIME-1111-1) :rate (s* (new-gname 'MMPRIME-0111-3) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0111-3) :to (new-gname 'LLPRIME-0111-2) :rate (s* (new-gname 'MMPRIME-0111-3) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0111-3) :to (new-gname 'LLPRIME-0111-3) :rate (s* (new-gname 'MMPRIME-0111-3) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0111-3) :to (new-gname 'LLPRIME-0111-4) :rate (s* (new-gname 'MMPRIME-0111-3) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-0111-3) :to (new-gname 'LLPRIME-1111-3) :rate (s* (new-gname 'LLPRIME-0111-3) (new-gname 'lam1) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-0111-4) :to (new-gname 'MMPRIME-1111-4) :rate (s* (new-gname 'MMPRIME-0111-4) (new-gname 'lam1) (new-gname 'zeta) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0111-4) :to (new-gname 'LLPRIME-1111-1) :rate (s* (new-gname 'MMPRIME-0111-4) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0111-4) :to (new-gname 'LLPRIME-0111-2) :rate (s* (new-gname 'MMPRIME-0111-4) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0111-4) :to (new-gname 'LLPRIME-0111-3) :rate (s* (new-gname 'MMPRIME-0111-4) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-0111-4) :to (new-gname 'LLPRIME-0111-4) :rate (s* (new-gname 'MMPRIME-0111-4) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'LLPRIME-0111-4) :to (new-gname 'LLPRIME-1111-4) :rate (s* (new-gname 'LLPRIME-0111-4) (new-gname 'lam1) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-1111-1) :to (new-gname 'LLPRIME-1111-1) :rate (s* (new-gname 'MMPRIME-1111-1) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1111-1) :to (new-gname 'LLPRIME-1111-2) :rate (s* (new-gname 'MMPRIME-1111-1) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1111-1) :to (new-gname 'LLPRIME-1111-3) :rate (s* (new-gname 'MMPRIME-1111-1) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1111-1) :to (new-gname 'LLPRIME-1111-4) :rate (s* (new-gname 'MMPRIME-1111-1) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-1111-2) :to (new-gname 'LLPRIME-1111-1) :rate (s* (new-gname 'MMPRIME-1111-2) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1111-2) :to (new-gname 'LLPRIME-1111-2) :rate (s* (new-gname 'MMPRIME-1111-2) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1111-2) :to (new-gname 'LLPRIME-1111-3) :rate (s* (new-gname 'MMPRIME-1111-2) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1111-2) :to (new-gname 'LLPRIME-1111-4) :rate (s* (new-gname 'MMPRIME-1111-2) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-1111-3) :to (new-gname 'LLPRIME-1111-1) :rate (s* (new-gname 'MMPRIME-1111-3) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1111-3) :to (new-gname 'LLPRIME-1111-2) :rate (s* (new-gname 'MMPRIME-1111-3) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1111-3) :to (new-gname 'LLPRIME-1111-3) :rate (s* (new-gname 'MMPRIME-1111-3) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1111-3) :to (new-gname 'LLPRIME-1111-4) :rate (s* (new-gname 'MMPRIME-1111-3) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'MMPRIME-1111-4) :to (new-gname 'LLPRIME-1111-1) :rate (s* (new-gname 'MMPRIME-1111-4) (new-gname 'lam1) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1111-4) :to (new-gname 'LLPRIME-1111-2) :rate (s* (new-gname 'MMPRIME-1111-4) (new-gname 'lam2) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1111-4) :to (new-gname 'LLPRIME-1111-3) :rate (s* (new-gname 'MMPRIME-1111-4) (new-gname 'lam3) (new-gname 'zeta) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'MMPRIME-1111-4) :to (new-gname 'LLPRIME-1111-4) :rate (s* (new-gname 'MMPRIME-1111-4) (new-gname 'lam4) (new-gname 'zeta) (new-gname 'pfast)))
@

RR states
<<>>=
(add-transition tbmodel   :from (new-gname 'RR-1000-1) :to (new-gname 'RR-1100-1) :rate (s* (new-gname 'RR-1000-1) (new-gname 'lam2) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'RR-1000-1) :to (new-gname 'RR-1010-1) :rate (s* (new-gname 'RR-1000-1) (new-gname 'lam3) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'RR-1000-1) :to (new-gname 'RR-1001-1) :rate (s* (new-gname 'RR-1000-1) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-0100-2) :to (new-gname 'RR-1100-2) :rate (s* (new-gname 'RR-0100-2) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'RR-0100-2) :to (new-gname 'RR-0110-2) :rate (s* (new-gname 'RR-0100-2) (new-gname 'lam3) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'RR-0100-2) :to (new-gname 'RR-0101-2) :rate (s* (new-gname 'RR-0100-2) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-0010-3) :to (new-gname 'RR-1010-3) :rate (s* (new-gname 'RR-0010-3) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'RR-0010-3) :to (new-gname 'RR-0110-3) :rate (s* (new-gname 'RR-0010-3) (new-gname 'lam2) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'RR-0010-3) :to (new-gname 'RR-0011-3) :rate (s* (new-gname 'RR-0010-3) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-0001-4) :to (new-gname 'RR-1001-4) :rate (s* (new-gname 'RR-0001-4) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'RR-0001-4) :to (new-gname 'RR-0101-4) :rate (s* (new-gname 'RR-0001-4) (new-gname 'lam2) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'RR-0001-4) :to (new-gname 'RR-0011-4) :rate (s* (new-gname 'RR-0001-4) (new-gname 'lam3) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-1100-1) :to (new-gname 'RR-1110-1) :rate (s* (new-gname 'RR-1100-1) (new-gname 'lam3) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'RR-1100-1) :to (new-gname 'RR-1101-1) :rate (s* (new-gname 'RR-1100-1) (new-gname 'lam4) (new-gname 'zeta) ))

(add-transition tbmodel   :from (new-gname 'RR-1100-2) :to (new-gname 'RR-1110-2) :rate (s* (new-gname 'RR-1100-2) (new-gname 'lam3) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'RR-1100-2) :to (new-gname 'RR-1101-2) :rate (s* (new-gname 'RR-1100-2) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-1010-1) :to (new-gname 'RR-1110-1) :rate (s* (new-gname 'RR-1010-1) (new-gname 'lam2) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'RR-1010-1) :to (new-gname 'RR-1011-1) :rate (s* (new-gname 'RR-1010-1) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-1010-3) :to (new-gname 'RR-1110-3) :rate (s* (new-gname 'RR-1010-3) (new-gname 'lam2) (new-gname 'zeta) ))
(add-transition tbmodel   :from (new-gname 'RR-1010-3) :to (new-gname 'RR-1011-3) :rate (s* (new-gname 'RR-1010-3) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-1001-1) :to (new-gname 'RR-1101-1) :rate (s* (new-gname 'RR-1001-1) (new-gname 'lam2) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'RR-1001-1) :to (new-gname 'RR-1011-1) :rate (s* (new-gname 'RR-1001-1) (new-gname 'lam3) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-1001-4) :to (new-gname 'RR-1101-4) :rate (s* (new-gname 'RR-1001-4) (new-gname 'lam2) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'RR-1001-4) :to (new-gname 'RR-1011-4) :rate (s* (new-gname 'RR-1001-4) (new-gname 'lam3) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-0110-2) :to (new-gname 'RR-1110-2) :rate (s* (new-gname 'RR-0110-2) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'RR-0110-2) :to (new-gname 'RR-0111-2) :rate (s* (new-gname 'RR-0110-2) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-0110-3) :to (new-gname 'RR-1110-3) :rate (s* (new-gname 'RR-0110-3) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'RR-0110-3) :to (new-gname 'RR-0111-3) :rate (s* (new-gname 'RR-0110-3) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-0101-2) :to (new-gname 'RR-1101-2) :rate (s* (new-gname 'RR-0101-2) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'RR-0101-2) :to (new-gname 'RR-0111-2) :rate (s* (new-gname 'RR-0101-2) (new-gname 'lam3) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-0101-4) :to (new-gname 'RR-1101-4) :rate (s* (new-gname 'RR-0101-4) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'RR-0101-4) :to (new-gname 'RR-0111-4) :rate (s* (new-gname 'RR-0101-4) (new-gname 'lam3) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-0011-3) :to (new-gname 'RR-1011-3) :rate (s* (new-gname 'RR-0011-3) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'RR-0011-3) :to (new-gname 'RR-0111-3) :rate (s* (new-gname 'RR-0011-3) (new-gname 'lam2) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-0011-4) :to (new-gname 'RR-1011-4) :rate (s* (new-gname 'RR-0011-4) (new-gname 'lam1) (new-gname 'zeta)))
(add-transition tbmodel   :from (new-gname 'RR-0011-4) :to (new-gname 'RR-0111-4) :rate (s* (new-gname 'RR-0011-4) (new-gname 'lam2) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-1110-1) :to (new-gname 'RR-1111-1) :rate (s* (new-gname 'RR-1110-1) (new-gname 'lam4) (new-gname 'zeta) ))

(add-transition tbmodel   :from (new-gname 'RR-1110-2) :to (new-gname 'RR-1111-2) :rate (s* (new-gname 'RR-1110-2) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-1110-3) :to (new-gname 'RR-1111-3) :rate (s* (new-gname 'RR-1110-3) (new-gname 'lam4) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-1101-1) :to (new-gname 'RR-1111-1) :rate (s* (new-gname 'RR-1101-1) (new-gname 'lam3) (new-gname 'zeta) ))

(add-transition tbmodel   :from (new-gname 'RR-1101-2) :to (new-gname 'RR-1111-2) :rate (s* (new-gname 'RR-1101-2) (new-gname 'lam3) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-1101-4) :to (new-gname 'RR-1111-4) :rate (s* (new-gname 'RR-1101-4) (new-gname 'lam3) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-1011-1) :to (new-gname 'RR-1111-1) :rate (s* (new-gname 'RR-1011-1) (new-gname 'lam2) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-1011-3) :to (new-gname 'RR-1111-3) :rate (s* (new-gname 'RR-1011-3) (new-gname 'lam2) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-1011-4) :to (new-gname 'RR-1111-4) :rate (s* (new-gname 'RR-1011-4) (new-gname 'lam2) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-0111-2) :to (new-gname 'RR-1111-2) :rate (s* (new-gname 'RR-0111-2) (new-gname 'lam1) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-0111-3) :to (new-gname 'RR-1111-3) :rate (s* (new-gname 'RR-0111-3) (new-gname 'lam1) (new-gname 'zeta)))

(add-transition tbmodel   :from (new-gname 'RR-0111-4) :to (new-gname 'RR-1111-4) :rate (s* (new-gname 'RR-0111-4) (new-gname 'lam1) (new-gname 'zeta)))
@

We do NOT currently have any transmission of resistant strains TO TB patients.

Infection the first time around:
one assumption is that only the predominant strain is transmitted.

<<>>=
;; Added SA 102912
(add-parameter tbmodel (new-gname 'zeta-uu))

(add-transition tbmodel   :from (new-gname 'UU) :to (new-gname 'MM-1000-1) :rate (s* (new-gname 'UU) (new-gname 'lam1)(new-gname 'zeta-uu) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'UU) :to (new-gname 'LL-1000-1) :rate (s* (new-gname 'UU) (new-gname 'lam1)(new-gname 'zeta-uu) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'ZZ) :to (new-gname 'MMPRIME-1000-1) :rate (s* (new-gname 'ZZ) (new-gname 'zeta) (new-gname 'lam1) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'ZZ) :to (new-gname 'LLPRIME-1000-1) :rate (s* (new-gname 'ZZ) (new-gname 'zeta) (new-gname 'lam1) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'QQ) :to (new-gname 'MM-1000-1) :rate (s* (new-gname 'QQ) (new-gname 'zeta) (new-gname 'lam1) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'QQ) :to (new-gname 'LL-1000-1) :rate (s* (new-gname 'QQ) (new-gname 'zeta) (new-gname 'lam1) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'UU) :to (new-gname 'MM-0100-2) :rate (s* (new-gname 'UU) (new-gname 'lam2)(new-gname 'zeta-uu) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'UU) :to (new-gname 'LL-0100-2) :rate (s* (new-gname 'UU) (new-gname 'lam2)(new-gname 'zeta-uu) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'ZZ) :to (new-gname 'MMPRIME-0100-2) :rate (s* (new-gname 'ZZ) (new-gname 'zeta) (new-gname 'lam2) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'ZZ) :to (new-gname 'LLPRIME-0100-2) :rate (s* (new-gname 'ZZ) (new-gname 'zeta) (new-gname 'lam2) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'QQ) :to (new-gname 'MM-0100-2) :rate (s* (new-gname 'QQ) (new-gname 'zeta) (new-gname 'lam2) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'QQ) :to (new-gname 'LL-0100-2) :rate (s* (new-gname 'QQ) (new-gname 'zeta) (new-gname 'lam2) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'UU) :to (new-gname 'MM-0010-3) :rate (s* (new-gname 'UU) (new-gname 'lam3)(new-gname 'zeta-uu) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'UU) :to (new-gname 'LL-0010-3) :rate (s* (new-gname 'UU) (new-gname 'lam3)(new-gname 'zeta-uu) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'ZZ) :to (new-gname 'MMPRIME-0010-3) :rate (s* (new-gname 'ZZ) (new-gname 'zeta) (new-gname 'lam3) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'ZZ) :to (new-gname 'LLPRIME-0010-3) :rate (s* (new-gname 'ZZ) (new-gname 'zeta) (new-gname 'lam3) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'QQ) :to (new-gname 'MM-0010-3) :rate (s* (new-gname 'QQ) (new-gname 'zeta) (new-gname 'lam3) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'QQ) :to (new-gname 'LL-0010-3) :rate (s* (new-gname 'QQ) (new-gname 'zeta) (new-gname 'lam3) (new-gname 'pfast)))

(add-transition tbmodel   :from (new-gname 'UU) :to (new-gname 'MM-0001-4) :rate (s* (new-gname 'UU) (new-gname 'lam4)(new-gname 'zeta-uu) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'UU) :to (new-gname 'LL-0001-4) :rate (s* (new-gname 'UU) (new-gname 'lam4)(new-gname 'zeta-uu) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'ZZ) :to (new-gname 'MMPRIME-0001-4) :rate (s* (new-gname 'ZZ) (new-gname 'zeta) (new-gname 'lam4) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'ZZ) :to (new-gname 'LLPRIME-0001-4) :rate (s* (new-gname 'ZZ) (new-gname 'zeta) (new-gname 'lam4) (new-gname 'pfast)))
(add-transition tbmodel   :from (new-gname 'QQ) :to (new-gname 'MM-0001-4) :rate (s* (new-gname 'QQ) (new-gname 'zeta) (new-gname 'lam4) (s- 1 (new-gname 'pfast))))
(add-transition tbmodel   :from (new-gname 'QQ) :to (new-gname 'LL-0001-4) :rate (s* (new-gname 'QQ) (new-gname 'zeta) (new-gname 'lam4) (new-gname 'pfast)))

(setf tbmodel-exogenous (deep-copy-model tbmodel))

(add-definition tbmodel (new-gname 'n-total) :value (apply #'s+ tb-state-list) :tag (new-gname 'n-total))
(add-parameter tbmodel (new-gname 'beta-smneg))
(add-parameter tbmodel (new-gname 'beta-smpos))
@

<<>>=
(add-definition tbmodel (new-gname 'smpos-1) :value (s+ (new-gname 'PP-1000-1) (new-gname 'PP-1100-1) (new-gname 'PP-1010-1) (new-gname 'PP-1001-1) (new-gname 'PP-1011-1) (new-gname 'PP-1101-1) (new-gname 'PP-1110-1) (new-gname 'PP-1111-1)(new-gname 'YYPOS-1000-1) (new-gname 'YYPOS-1100-1) (new-gname 'YYPOS-1010-1) (new-gname 'YYPOS-1001-1) (new-gname 'YYPOS-1011-1) (new-gname 'YYPOS-1101-1) (new-gname 'YYPOS-1110-1) (new-gname 'YYPOS-1111-1)(new-gname 'WWPOS-1000-1) (new-gname 'WWPOS-1100-1) (new-gname 'WWPOS-1010-1) (new-gname 'WWPOS-1001-1) (new-gname 'WWPOS-1011-1) (new-gname 'WWPOS-1101-1) (new-gname 'WWPOS-1110-1) (new-gname 'WWPOS-1111-1)) :tag (new-gname 'smearpositive1))
(add-definition tbmodel (new-gname 'smpos-2) :value (s+ (new-gname 'PP-0100-2) (new-gname 'PP-1100-2) (new-gname 'PP-0110-2) (new-gname 'PP-0101-2) (new-gname 'PP-0111-2) (new-gname 'PP-1101-2) (new-gname 'PP-1110-2) (new-gname 'PP-1111-2)(new-gname 'YYPOS-0100-2) (new-gname 'YYPOS-1100-2) (new-gname 'YYPOS-0110-2) (new-gname 'YYPOS-0101-2) (new-gname 'YYPOS-0111-2) (new-gname 'YYPOS-1101-2) (new-gname 'YYPOS-1110-2) (new-gname 'YYPOS-1111-2)(new-gname 'WWPOS-0100-2) (new-gname 'WWPOS-1100-2) (new-gname 'WWPOS-0110-2) (new-gname 'WWPOS-0101-2) (new-gname 'WWPOS-0111-2) (new-gname 'WWPOS-1101-2) (new-gname 'WWPOS-1110-2) (new-gname 'WWPOS-1111-2)) :tag (new-gname 'smearpositive2))
(add-definition tbmodel (new-gname 'smpos-3) :value (s+ (new-gname 'PP-0010-3) (new-gname 'PP-1010-3) (new-gname 'PP-0110-3) (new-gname 'PP-0011-3) (new-gname 'PP-0111-3) (new-gname 'PP-1011-3) (new-gname 'PP-1110-3) (new-gname 'PP-1111-3)(new-gname 'YYPOS-0010-3) (new-gname 'YYPOS-1010-3) (new-gname 'YYPOS-0110-3) (new-gname 'YYPOS-0011-3) (new-gname 'YYPOS-0111-3) (new-gname 'YYPOS-1011-3) (new-gname 'YYPOS-1110-3) (new-gname 'YYPOS-1111-3)(new-gname 'WWPOS-0010-3) (new-gname 'WWPOS-1010-3) (new-gname 'WWPOS-0110-3) (new-gname 'WWPOS-0011-3) (new-gname 'WWPOS-0111-3) (new-gname 'WWPOS-1011-3) (new-gname 'WWPOS-1110-3) (new-gname 'WWPOS-1111-3)) :tag (new-gname 'smearpositive3))
(add-definition tbmodel (new-gname 'smpos-4) :value (s+ (new-gname 'PP-0001-4) (new-gname 'PP-1001-4) (new-gname 'PP-0101-4) (new-gname 'PP-0011-4) (new-gname 'PP-0111-4) (new-gname 'PP-1011-4) (new-gname 'PP-1101-4) (new-gname 'PP-1111-4)(new-gname 'YYPOS-0001-4) (new-gname 'YYPOS-1001-4) (new-gname 'YYPOS-0101-4) (new-gname 'YYPOS-0011-4) (new-gname 'YYPOS-0111-4) (new-gname 'YYPOS-1011-4) (new-gname 'YYPOS-1101-4) (new-gname 'YYPOS-1111-4)(new-gname 'WWPOS-0001-4) (new-gname 'WWPOS-1001-4) (new-gname 'WWPOS-0101-4) (new-gname 'WWPOS-0011-4) (new-gname 'WWPOS-0111-4) (new-gname 'WWPOS-1011-4) (new-gname 'WWPOS-1101-4) (new-gname 'WWPOS-1111-4)) :tag (new-gname 'smearpositive4))
(add-definition tbmodel (new-gname 'smneg-1) :value (s+ (new-gname 'SS-1000-1) (new-gname 'SS-1100-1) (new-gname 'SS-1010-1) (new-gname 'SS-1001-1) (new-gname 'SS-1011-1) (new-gname 'SS-1101-1) (new-gname 'SS-1110-1) (new-gname 'SS-1111-1)) :tag (new-gname 'smearnegative1))
(add-definition tbmodel (new-gname 'smneg-2) :value (s+ (new-gname 'SS-0100-2) (new-gname 'SS-1100-2) (new-gname 'SS-0110-2) (new-gname 'SS-0101-2) (new-gname 'SS-0111-2) (new-gname 'SS-1101-2) (new-gname 'SS-1110-2) (new-gname 'SS-1111-2)) :tag (new-gname 'smearnegative2))
(add-definition tbmodel (new-gname 'smneg-3) :value (s+ (new-gname 'SS-0010-3) (new-gname 'SS-1010-3) (new-gname 'SS-0110-3) (new-gname 'SS-0011-3) (new-gname 'SS-0111-3) (new-gname 'SS-1011-3) (new-gname 'SS-1110-3) (new-gname 'SS-1111-3)) :tag (new-gname 'smearnegative3))
(add-definition tbmodel (new-gname 'smneg-4) :value (s+ (new-gname 'SS-0001-4) (new-gname 'SS-1001-4) (new-gname 'SS-0101-4) (new-gname 'SS-0011-4) (new-gname 'SS-0111-4) (new-gname 'SS-1011-4) (new-gname 'SS-1101-4) (new-gname 'SS-1111-4)) :tag (new-gname 'smearnegative4))
@

<<>>=
(redefine tbmodel (new-gname 'lam1) :value (s+ (s/ (s* (new-gname 'beta-smpos) (new-gname 'smpos-1)) (new-gname 'n-total)) (s/ (s* (new-gname 'beta-smneg) (new-gname 'smneg-1)) (new-gname 'n-total))))
(redefine tbmodel (new-gname 'lam2) :value (s+ (s/ (s* (new-gname 'beta-smpos) (new-gname 'smpos-2)) (new-gname 'n-total)) (s/ (s* (new-gname 'beta-smneg) (new-gname 'smneg-2)) (new-gname 'n-total))))
(redefine tbmodel (new-gname 'lam3) :value (s+ (s/ (s* (new-gname 'beta-smpos) (new-gname 'smpos-3)) (new-gname 'n-total)) (s/ (s* (new-gname 'beta-smneg) (new-gname 'smneg-3)) (new-gname 'n-total))))
(redefine tbmodel (new-gname 'lam4) :value (s+ (s/ (s* (new-gname 'beta-smpos) (new-gname 'smpos-4)) (new-gname 'n-total)) (s/ (s* (new-gname 'beta-smneg) (new-gname 'smneg-4)) (new-gname 'n-total))))
@

<<>>=
(add-state tbmodel (new-gname 'cumul-smneg-incidence))
(add-state tbmodel (new-gname 'cumul-smpos-incidence))
@

<<>>=
; collect up all transitions which have to do with smpos incidence
;   exclude transitions from smneg to smpos or treatment discontinuation, etc.
;   new-fast-latent-states  
;   new-slow-latent-states  
;   smearnegative-states
;   smearpositive-states
;   old-fast-latent-states
;   old-slow-latent-states
;   incomplete-cure-states
;   treated-smneg-states-standard
;   treated-smneg-states-specialized
;   treated-smpos-states-standard
;   treated-smpos-states-specialized

(defvar notnew (append treated-smneg-states-standard treated-smneg-states-specialized treated-smpos-states-standard treated-smpos-states-specialized smearnegative-states smearpositive-states))

(defun smpos-p (ss) (member ss smearpositive-states))
(defun smneg-p (ss) (member ss smearnegative-states))
(defun incid-smpos-p (tt)
  (and (smpos-p (transition-to tt))
       (not (member (transition-from tt) notnew))))
(defun incid-smneg-p (tt)
  (and (smneg-p (transition-to tt))
       (not (member (transition-from tt) notnew))))
@

<<>>=
(dolist (tt (model-transitions tbmodel))
  (if (incid-smpos-p tt)
      (add-transition tbmodel  :from (new-gname '_SOURCE) :to (new-gname 'cumul-smpos-incidence) :rate (transition-value tt))))
(dolist (tt (model-transitions tbmodel))
  (if (incid-smneg-p tt)
      (add-transition tbmodel  :from (new-gname '_SOURCE) :to (new-gname 'cumul-smneg-incidence) :rate (transition-value tt))))
@

<<>>=
; age model
(defvar agemodel (new-model 'agemodel))
(add-flow-series agemodel :states (indexed (new-gname 'AGE) 0 19) :parameter (new-gname 'agingrate) :stratify NIL)
(let ((tmpname NIL))
  (dolist (ss (all-keys (model-states agemodel)))
    (setf tmpname (merge-names (new-gname 'mu) ss))
    (add-parameter agemodel tmpname)
    (add-transition agemodel  :from ss :to (new-gname '_SINK) :rate (s* ss tmpname))))
@

<<>>=
(defvar littleagemodel (new-model 'littleagemodel))
(add-flow-series littleagemodel :states (indexed (new-gname 'AGE) 0 4) :parameter (new-gname 'agingrate) :stratify T)
(add-parameter littleagemodel (new-gname 'biglam))
(add-transition littleagemodel :from (new-gname '_SOURCE) :to (@ (new-gname 'AGE) 0) :rate (new-gname 'biglam) )
@

<<>>=
(let ((tmpname NIL))
  (dolist (ss (all-keys (model-states littleagemodel)))
    (setf tmpname (merge-names (new-gname 'mu) ss))
    (add-parameter littleagemodel tmpname)
    (add-transition littleagemodel  :from ss :to (new-gname '_SINK) :rate (s* ss tmpname))))
@

<<>>=
; nutrition model
(defvar nutrimodel (new-model 'nutrimodel))
(add-state nutrimodel (new-gname 'MALNUT))
(add-state nutrimodel (new-gname 'NOUR))
(add-parameter nutrimodel (new-gname 'starve))
(add-parameter nutrimodel (new-gname 'feed))
@

<<>>=
; immigration model
(defvar migmodel (new-model 'migmodel))
(add-state migmodel (new-gname 'USAB))
(add-state migmodel (new-gname 'FB))
@

<<>>=
(defvar hivprogmodel (new-model 'hivprogmodel))
@

<<>>=
(setf dictry (new-dictionary))
(suppress-dependence dictry :of (new-gname 'starve))
(suppress-dependence dictry :of (new-gname 'feed))
(suppress-dependence dictry :of (new-gname 'lam1))
(suppress-dependence dictry :of (new-gname 'lam2))
(suppress-dependence dictry :of (new-gname 'lam3))
(suppress-dependence dictry :of (new-gname 'lam4))
(suppress-dependence dictry :of (new-gname 'beta-smpos))
(suppress-dependence dictry :of (new-gname 'beta-smneg))
(defvar tbnmodel (cross tbmodel nutrimodel dictry))
@

<<>>=
(setf agetbd (new-dictionary))
(suppress-dependence agetbd :of (new-gname 'agingrate)) 
(suppress-dependence agetbd :of (pa agemodel))
(suppress-dependence agetbd :of (remove 'smpos-prob (pa tbmodel)))
;(suppress-dependence agetbd :of (new-gname 'n-total))
(suppress-dependence agetbd :of (new-gname 'lam1))
(suppress-dependence agetbd :of (new-gname 'lam2))
(suppress-dependence agetbd :of (new-gname 'lam3))
(suppress-dependence agetbd :of (new-gname 'lam4))
(suppress-dependence dictry :of (new-gname 'beta-smpos))
(suppress-dependence dictry :of (new-gname 'beta-smneg))
;(defvar model2 (cross tbmodel agemodel agetbd))
@

<<>>=
(defvar nkmodel (cross tbnmodel littleagemodel agetbd))
@

Now compile:
<<>>=
(finish nkmodel)
@

\vfill
\end{document}
