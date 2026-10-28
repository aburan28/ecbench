"""Artin-Schreier trace-block elimination experiment (rank 4). Reuses certified harness; modifies nothing there."""
import sys, os, random, itertools, subprocess, tempfile, re, hashlib, time, json
HERE = os.path.dirname(os.path.abspath(__file__))
HARN = os.path.abspath(os.path.join(HERE, '..', '..', 'harness'))
sys.path.insert(0, HARN)
import gf2, boolean_weil as bw, pdp, exp_chain as EC
from gf2 import gmul, gsq, bin_add, bin_points

A2, A6 = 1, 1
PRE = ['nice', '-n', '10', 'timeout', '300', 'prlimit', '--as=2500000000']

# ---------- field helpers ----------
def tr(x, n):
    t, y = 0, x
    for _ in range(n):
        t ^= y; y = gsq(y, n)
    assert t in (0, 1)
    return t

def coords(basis, val):
    """coordinates of val in F2-basis (list of ints); None if not in span."""
    piv = {}
    for i, b in enumerate(basis):
        v, mk = b, 1 << i
        while v:
            h = v.bit_length() - 1
            if h in piv: v ^= piv[h][0]; mk ^= piv[h][1]
            else: piv[h] = (v, mk); break
    mk = 0; v = val
    while v:
        h = v.bit_length() - 1
        if h not in piv: return None
        v ^= piv[h][0]; mk ^= piv[h][1]
    return [(mk >> i) & 1 for i in range(len(basis))]

def span(basis):
    S = [0]
    for b in basis: S += [s ^ b for s in S]
    return S

def adapt_V(B, n):
    trs = [tr(b, n) for b in B]
    if not any(trs): return None
    j = trs.index(1); b1 = B[j]
    return [b1] + [B[i] ^ (b1 if trs[i] else 0) for i in range(len(B)) if i != j]

def adapt_E(n):
    j0 = next(k for k in range(n) if tr(1 << k, n))
    e0 = 1 << j0
    return [e0] + [(1 << k) ^ (e0 if tr(1 << k, n) else 0) for k in range(n) if k != j0]

# ---------- numeric systems ----------
def S3n(x1, x2, x3, n): return gf2.S3_bin(x1, x2, x3, A6, n)
def S4n(x1, x2, x3, x4, n):
    def co(a, b):
        B = gmul(a, b, n); return gsq(a ^ b, n), B, gsq(B, n) ^ A6
    a1, b1, c1 = co(x1, x2); a2, b2, c2 = co(x3, x4)
    t1 = gsq(gmul(a1, c2, n) ^ gmul(a2, c1, n), n)
    t2 = gmul(gmul(a1, b2, n) ^ gmul(a2, b1, n), gmul(b1, c2, n) ^ gmul(b2, c1, n), n)
    return t1 ^ t2

# ---------- symbolic helpers ----------
def subst(eqs, assign):
    out = []
    for f in eqs:
        g = set()
        for mon in f:
            if any(v in assign and assign[v] == 0 for v in mon): continue
            m2 = frozenset(v for v in mon if v not in assign)
            if m2 in g: g.discard(m2)
            else: g.add(m2)
        if g: out.append(g)
    return out

def evalpoly(f, asg):
    s = 0
    for mon in f:
        p = 1
        for v in mon: p &= asg[v]
        s ^= p
    return s

def maxdeg_in(eqs): return max((len(m) for f in eqs for m in f), default=0)

# ---------- cell construction ----------
def make_cell(n, l, m, model, seed0=1):
    """returns dict with everything; bumps seed until trace block non-vacuous and decomposition exists."""
    seed = seed0
    while True:
        rng = random.Random(seed)
        B = pdp.random_subspace_basis(n, l, rng)
        dec = pdp.find_decomposition(A2, A6, n, B, m, rng)
        Bad = adapt_V(B, n)
        if dec is not None and Bad is not None and dec['xR'] != 0: break
        seed += 1
    xR = dec['xR']; xs = dec['xs']
    Eb = adapt_E(n)
    split = (model == 'split')
    def build(Bx, Ex):
        X = []; allv = []
        for i in range(m):
            xi, vi = bw.var_elt('x', i + 1, Bx, n); X.append(xi); allv += vi
        XR = EC.const_elt(xR, n)
        if m == 2:
            eqs = bw.weil_equations(bw.S3_bin_sym(X[0], X[1], XR, A6, n), n)
        elif model == 'unsplit':
            eqs = bw.weil_equations(bw.S4_unsplit(X[0], X[1], X[2], XR, A6, n), n)
        else:
            if Ex is None: Aa, av = EC.free_elt('a', 1, n)
            else: Aa, av = bw.var_elt('a', 1, Ex, n)
            allv += av
            eqs = bw.weil_equations(bw.S3_bin_sym(X[0], X[1], Aa, A6, n), n) + \
                  bw.weil_equations(bw.S3_bin_sym(Aa, X[2], XR, A6, n), n)
        return eqs, sorted(set(allv))
    eqs, vars_ = build(Bad, Eb if split else None)
    eqs0, vars0 = build(B, None)      # original harness coordinates (reference)
    # planted assignment in adapted coordinates
    asg = {}
    for i, x in enumerate(xs):
        c = coords(Bad, x)
        for j, b in enumerate(c): asg[('x', i + 1, j)] = b
    aux_val = None
    if split:
        pts = [P for P in bin_points(A2, A6, n) if P is not None]
        # recover planted points to get x(P1+P2)
        P = {p[0]: p for p in pts}
        # find y's: re-derive via decomposition list
        # find_decomposition returned only xs; locate points by brute force over sign choices
        for sg in itertools.product(*[[p for p in pts if p[0] == x] for x in xs]):
            S = None
            for q in sg: S = bin_add(S, q, A2, A6, n)
            if S is not None and S[0] == xR:
                S12 = bin_add(sg[0], sg[1], A2, A6, n)
                aux_val = S12[0]; break
        assert aux_val is not None
        c = coords(Eb, aux_val)
        for k, b in enumerate(c): asg[('a', 1, k)] = b
    for f in eqs: assert evalpoly(f, asg) == 0, "planted not a solution"
    return dict(n=n, l=l, m=m, model=model, seed=seed, B=B, Bad=Bad, Eb=Eb, xs=xs, xR=xR, aux=aux_val,
                eqs=eqs, vars=vars_, eqs0=eqs0, vars0=vars0, planted=asg, split=split)

# ---------- numeric solution sets / brute force ----------
def brute(cell):
    n, l, m = cell['n'], cell['l'], cell['m']; xR = cell['xR']; Bad = cell['Bad']
    V = span(Bad)
    pts = [P for P in bin_points(A2, A6, n) if P is not None]
    fb = [P for P in pts if P[0] in set(V)]
    BF = set(); BFc = set()
    for tup in itertools.product(fb, repeat=m):
        S = None
        for q in tup: S = bin_add(S, q, A2, A6, n)
        if S is not None and S[0] == xR:
            BF.add(tuple(q[0] for q in tup))
            if m == 2 or bin_add(tup[0], tup[1], A2, A6, n) is not None: BFc.add(tuple(q[0] for q in tup))
    sols = []   # list of (xtuple, auxval or None)
    if m == 2:
        for x1 in V:
            for x2 in V:
                if S3n(x1, x2, xR, n) == 0: sols.append(((x1, x2), None))
    elif cell['model'] == 'unsplit':
        for t in itertools.product(V, repeat=3):
            if S4n(t[0], t[1], t[2], xR, n) == 0: sols.append((t, None))
    else:
        for x1 in V:
            for x2 in V:
                As = [a for a in range(1 << n) if S3n(x1, x2, a, n) == 0]
                for x3 in V:
                    for a in As:
                        if S3n(a, x3, xR, n) == 0: sols.append(((x1, x2, x3), a))
    SFx = set(s[0] for s in sols)
    return dict(BF=BF, BFc=BFc, sols=sols, SFx=SFx)

def branch_key(cell, sol, arm):
    xt, a = sol; n, l = cell['n'], cell['l']; bits = []
    for x in xt:
        bits.append(tr(x, n) if arm == 'as' else coords(cell['Bad'], x)[l - 1])
    if cell['split']:
        bits.append(tr(a, n) if arm == 'as' else coords(cell['Eb'], a)[n - 1])
    return tuple(bits)

def branch_assigns(cell, arm):
    n, l, m = cell['n'], cell['l'], cell['m']
    names = [('x', i + 1, 0 if arm == 'as' else l - 1) for i in range(m)]
    if cell['split']: names.append(('a', 1, 0 if arm == 'as' else n - 1))
    for bits in itertools.product([0, 1], repeat=len(names)):
        yield bits, dict(zip(names, bits))

# ---------- engines (nice/timeout/prlimit) ----------
def singular(eqs, vars_):
    vs = sorted(vars_); vmap = {v: f"c{i}" for i, v in enumerate(vs)}
    gens = [pdp.to_singular_poly(f, vmap) for f in eqs] + [f"{vmap[v]}^2+{vmap[v]}" for v in vs]
    script = (f"ring r=2,({','.join(vmap[v] for v in vs)}),dp;\noption(prot);\nideal I={','.join(gens)};\n"
              "ideal G=std(I);\nprint(\"DIM \"+string(dim(G)));\nprint(\"VDIM \"+string(vdim(G)));\n$\n")
    t0 = time.time()
    try:
        p = subprocess.run(PRE + ['Singular', '-q'], input=script, capture_output=True, text=True, timeout=330)
    except subprocess.TimeoutExpired:
        return dict(status='timeout', deg=None, vdim=None, dim=None)
    out = p.stdout + p.stderr
    if p.returncode != 0 or 'VDIM' not in out:
        return dict(status='timeout_or_fail rc=%d' % p.returncode, deg=None, vdim=None, dim=None)
    r = pdp.parse_singular_degree(out)
    return dict(status='ok', deg=r['solving_degree'], vdim=r['vdim'], dim=r['dim'])

def msolve(eqs, vars_):
    inp = __import__('exp_engine_compare').to_msolve(eqs, vars_)
    with tempfile.NamedTemporaryFile('w', suffix='.ms', delete=False) as fh:
        fh.write(inp); path = fh.name
    try:
        try:
            p = subprocess.run(PRE + ['msolve', '-t', '1', '-v', '2', '-f', path, '-g', '2', '-o', '/dev/null'],
                               capture_output=True, text=True, timeout=330)
        except subprocess.TimeoutExpired:
            return dict(status='timeout', deg=None)
        out = p.stdout + p.stderr
        degs = [int(m.group(1)) for line in out.splitlines()
                for m in [re.match(r'\s*(\d+)\s+\d+\s+\d+\s+\d+ x \d+', line)] if m]
        if p.returncode != 0: return dict(status='fail rc=%d' % p.returncode, deg=None)
        return dict(status='ok', deg=max(degs) if degs else None)
    finally:
        os.unlink(path)

def macaulay(eqs, vars_):
    from math import comb
    v = len(vars_); D = 0
    for d in range(1, 6):
        if sum(comb(v, k) for k in range(d + 1)) <= 4500: D = d
    if D < 3: return dict(ff=None, Dmax=D, note='skipped')
    mac = pdp.macaulay_rank_seq(eqs, vars_, D)
    return dict(ff=mac['first_fall_degree'], Dmax=D)

def sysh(eqs, vars_): return EC.sys_hash(eqs, vars_)

def measure(eqs, vars_, mac=False):
    r = dict(nvars=len(vars_), neqs=len(eqs), din=maxdeg_in(eqs), sha=sysh(eqs, vars_))
    r['singular'] = singular(eqs, vars_); r['msolve'] = msolve(eqs, vars_)
    if mac: r['mac'] = macaulay(eqs, vars_)
    return r
