# PROTOCOL (written before running; 2026-10-23) -- rank 5: degree-2 / Kummer-style encoding

## Hypothesis
Replacing each S_3 link by its Vieta ("sum and difference abscissa") form
  D*(s+t) = 2((AB+a)(A+B)+2b),  D*s*t = (AB-a)^2-4b(A+B),  D=(A-B)^2,
where s=x(P+Q), t=x(P-Q), A=x(P), B=x(Q), and optionally flattening further to a fully
quadratic system with symmetric-function auxiliaries (e=A+B, w=AB, delta=e^2-4w, sigma=s+t, pi=st),
lowers the SOLVING degree (max degree reached in F4 / std) below the S_3-chain. Prior: it is neutral or worse.

## Setting
Prime field F_p (p<=2017), y^2=x^3+ax+b, factor base V={P: x(P)^d=c}, c a d-th power. m=3 main, m=2 control.
Target x_R = x(P1+P2+P3) for a planted decomposition with x(Pi) in V. All systems include x_i^d=c (i<=m).

## Systems compared (m=3; unknowns x1,x2,x3 + auxiliaries)
(a) unsplit: S4(x1,x2,x3,xR)=0, S4=Res_u(S3(x1,x2,u),S3(u,x3,xR)).
(b) chain: S3(x1,x2,u)=0, S3(u,x3,xR)=0.
(c1) Vieta pairs: link(x1,x2;u,u'), link(u,x3;xR,v)  (2 Vieta eqs per link; u',v are the "difference" auxiliaries).
(c2) fully quadratic flattening of c1 (e,w,delta,sigma,pi auxiliaries per link, s,t from sigma,pi).
m=2: (a)=(b)=S3(x1,x2,xR); (c1) link(x1,x2;xR,v); (c2) its flattening.

## Correctness (before degrees)
Brute force on each curve: for all (x1,x2,x3) in V^3 with all xi on the curve (f(x) square), compare
the F_p-solvability of each system's projection to (x1,x2,x3) against ground truth
"exists signs: x(+-P1+-P2+-P3)=xR" computed by the group law. Report false pos/neg (must be 0).
Planted decomposition: ideal + (x_i - x_i*) must not be <1> (Singular), for every system.

## Measurements
msolve F4 (-t 1 -v 2, max degree column), Singular std (dp, option(prot), max [k:d]), vdim (Singular),
#vars, #eqs, max total degree, Bezout = product of the n smallest generator degrees (n=#vars; the
Bezout bound for overdetermined systems; also the full product of all gens is stored). Run twice, diff hashes/degrees/vdim.

## Falsifier (fixed; from survey_E rank 5, sharpened)
Lever is FALSIFIED if for (c1) and (c2) [each], d_solv(c) >= d_solv(b) on BOTH engines with no Bezout drop
(Bezout(c) >= Bezout(b)) -- on the majority of configs. "Lowers" requires d_solv(c) < d_solv(b) on both engines.
Mixed engines => inconclusive/engine-dependent. Outcomes: lowers / neutral (= equal) / raises.
Caveat declared in advance: total degree of Vieta link is 4 (same as S3) in c1; only c2 is genuinely degree 2/3 per equation.

## Amendment 2026-10-23 (after first pass, falsifier untouched)
First pass parsed Singular `[k:d]` tokens, which only encode the input ([65535:2]); inspection of the protocol
showed the step degrees are the bare integers (4,6,...,23) between the s/- markers ("(n)" = pair-set sizes). The
Singular degree is therefore the max bare integer in the protocol. First-pass Singular numbers (2..5) are discarded.
Also: planted x1,x2,x3 are a decomposition only for m=3; for m=2 the "planted_vdim=0" is expected and not a check.

## Amendment 2 (2026-10-23, per coordinator; falsifier unchanged)
Singular d_solv := least D such that std(I) under degBound=D already has vdim == vdim of the unbounded std(I)
(identical to rank3 run_solvers.py; scan D=1..60, 300 s cap -> missing cell recorded as "not reached").
The previous bare-integer parse is discarded. Planted check = (G + (x_i - x_i*)) != <1> (m=3; for m=2 the target is
re-derived as x(P1+P2) from the same two planted points). Reported for each system: raw d_solv and offset d_solv - d.
Added factor-base variants d in {4,6} (p=1009) to expose the offset; configs: (1009,4),(1009,6),(1009,8),(1009,16),(2017,8).
Bezout reported both as product of the n smallest generator degrees and as full product.
