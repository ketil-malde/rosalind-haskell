Proposed Haskell solutions to Rosalind problems
===============================================

Recently, I attended the "Bioinformatics in Bergen" meeting. Part of
this was a small (two hours) hackathon/coding contest, where we were
to solve various problems from the
[Rosalind](http://rosalind.info/problems/list-view/) set.

I decided to solve as many as possible using only Haskell, no
libraries allowed outside of the Prelude.  (Admittedly, I used some
list functions from the standard library, but I later implemented them
to eliminate this dependency, and they are now at the bottom of the
source file.)

The main function reads test data from the
[BIB 2016](https://github.com/PierreBedoucha/bib2016) repository and produces
the (or so I hope!) correct output.

I hope this can be of use to people who know bioinformatics, and are
interested in Haskell - or vice versa.  Needless to say, the
implementations are naïve and probably scale less than ideal.

