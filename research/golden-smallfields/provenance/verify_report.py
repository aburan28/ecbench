"""Independent consistency checks for the report's numerical records."""
import json
import math
from pathlib import Path
from sympy import factorint, isprime

ROOT = Path(__file__).resolve().parent

def factors_value(factors):
    return math.prod(int(p) ** int(e) for p, e in factors.items())

def verify_certificate(cert):
    n = int(cert['n'])
    if n == 2:
        return True
    f = cert.get('n_minus_1_factors', cert.get('factors'))
    if not f:
        return bool(isprime(n))
    assert factors_value(f) == n - 1
    children = cert.get('prime_factor_certificates', {})
    for q in f:
        q = int(q)
        if str(q) in children:
            verify_certificate(children[str(q)])
        else:
            assert isprime(q)
    witnesses = cert.get('witness_for_each_prime_factor', cert.get('witnesses', {}))
    for q in f:
        q = int(q)
        assert str(q) in witnesses
        a = int(witnesses[str(q)])
        assert pow(a, n-1, n) == 1
        assert math.gcd(pow(a, (n-1)//q, n)-1, n) == 1
    return True

def verify_curve(c):
    q, t, n, r, h, k = (int(c[x]) for x in ['q','t','N','r','h','k'])
    assert n == q + 1 - t
    assert t*t <= 4*q
    assert int(c['Dpi']) == t*t - 4*q
    assert factors_value(c['factors']) == n
    assert all(isprime(int(p)) for p in c['factors'])
    assert r == max(int(p) for p in c['factors'])
    assert r*h == n
    assert math.gcd(q,r) == 1
    assert pow(q,k,r) == 1
    for ell in factorint(k):
        assert pow(q,k//int(ell),r) != 1
    d1,d2 = map(int,c['group'])
    assert d1*d2 == n and d2%d1 == 0 and (q-1)%d1 == 0
    if c.get('DK') is not None and c.get('fpi') is not None:
        assert int(c['Dpi']) == int(c['DK'])*int(c['fpi'])**2
    if c.get('ordinary') and c.get('Dend') is not None:
        ratio = int(c['Dend'])//int(c['DK'])
        assert int(c['Dend']) == int(c['DK'])*math.isqrt(ratio)**2
        assert int(c['fpi'])%math.isqrt(ratio) == 0

def main():
    data = json.loads((ROOT/'normalized_data.json').read_text())
    for c in data['curves']:
        verify_curve(c)
    lookup = {c['id']:c for c in data['curves']}
    for c in data['curves']:
        if c.get('twist_id') in lookup:
            d = lookup[c['twist_id']]
            assert c['q'] == d['q'] and c['t'] == -d['t']
            assert c['N'] + d['N'] == 2*(c['q']+1)
            assert c['Dpi'] == d['Dpi']
    certificates = json.loads((ROOT/'prime_certificates.json').read_text())
    for value, cert in certificates.items():
        n=int(value)
        assert n == cert['n']
        if n == 2:
            continue
        f=cert['n_minus_1_factors']
        assert factors_value(f) == n-1
        a=cert['a']
        assert pow(a,n-1,n) == 1
        for prime in f:
            assert str(prime) in certificates
            assert math.gcd(pow(a,(n-1)//int(prime),n)-1,n) == 1
    binary_certificates=json.loads((ROOT/'binary_certificates.json').read_text())
    for cert in binary_certificates.values():
        verify_certificate(cert)
    print(f"Validated {len(certificates)} flat prime-field Lucas certificates and {len(binary_certificates)} recursive binary certificates.")
    print(f"Passed independent identities, factor primality, group structure, and exact embedding-degree checks for {len(lookup)} curves.")

if __name__ == '__main__':
    main()
