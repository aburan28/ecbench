# PROTOCOL — rank 6: Lattès / isogeny layer maps in a PKM-style tower factor base (written before running)

Academic toy measurement of Gröbner solving degree. No claims about deployed curves.

## Hypothesis
PKM factor base F = {x in F_p : L(x) in C}, L = L_k o ... o L_1, |C| = s. Choosing the layers L_j with
arithmetic structure (Lattès [2]-map of an auxiliary curve, degree 4; x-maps of rational 2-/3-/(cyclic 4)-isogenies,
via Velu) might lower the degree reached by Groebner engines on the decomposition system relative to
generic rational maps of the same degree. Survey_E predicts probably not (decomposition counts ~ random).

## Systems (per prime p, target curve E: y^2=x^3+ax+b, m=2)
Unknowns y_{i,j}, i in {1,2}, j in 0..k. Equations:
 - layer: y_{i,j+1}*den_j(y_{i,j}) - num_j(y_{i,j}) = 0  (deg d+1)
 - base:  prod_{c in C}(y_{i,k}-c) = 0  (deg s)
 - S_3(y_{1,0}, y_{2,0}, xR) = 0, xR a constant (the target's x-coordinate), planted so that x(P1),x(P2) in F.
Families at EQUAL (k, d, s, #vars): lattes (d=4,k=2; L_j = x([2]) of one aux curve E', so L = x([4])),
isogeny (Velu chains of rational isogenies: d=2,k=3; d=3,k=2; cyclic-4 = two 2-isogenies per layer, d=4,k=2),
generic (random coprime num/den, degree d per layer), kummer (x -> x^d per layer). Control: "random base"
= random x-set of the same size as the tower base (no tower, only counts decompositions).
C = s random elements of the image L(F_p) (so every c has non-empty fibre), s ~ sqrt(p).
Seeds fixed; the system hash is recorded.

## Measurements
Per system: variables, Bezout (prod of degrees), msolve F4 max degree (-t 1 -v 2 round table), Singular std (dp)
max [k:d] degree and vdim, factor-base size |F|, average number of ordered decompositions of a random target
(over 100 targets with points) vs the random-base control, planted decomposition recovered (Singular: planted
point not eliminated, plus direct root check of all equations). Correctness: S_3 vs brute-force on p=101
(0 false pos / neg) and isogeny maps checked to send curve points to curve points.
Each deterministic measurement is run twice (run1/run2) and diffed. Solver wrapper:
nice -n 10 timeout 300 prlimit --as=2500000000. Timeout = recorded as ">= d, not reached".

## Falsifier (fixed before running; from survey_E §5 rank 6)
For a configuration (p,k,d), a structured family (lattes or isogeny) is "not better" if its d_solv equals the
generic family's d_solv. Falsified (verdict NEUTRAL) if, on BOTH engines, structured d_solv == generic d_solv
(allowing +-0) over all configurations measured. LOWERS only if structured < generic on both engines in a
majority of configurations with generic-family spread (over 3 seeds) not covering the difference.
RAISES if structured > generic likewise. Otherwise INCONCLUSIVE. Generic = median over 3 seeds.

## Not done up front
m=3 (S_3 chain) only if time permits; otherwise reported tier T5.
