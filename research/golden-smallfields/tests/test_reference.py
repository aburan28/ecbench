"""Meaningful toy-curve conformance and release tamper tests."""
import copy
import unittest
from ecbench.reference import Field,Curve,ECError
from ecbench.validate import load_dataset,verify,certify


def toy(kind):
    binary=kind=='binary'
    f=Field({'kind':kind,'cardinality':str(32 if binary else 17),'element_bits':5,'coordinate_bytes':1,
             'representation':{'modulus':'0x25' if binary else '0x11'}})
    values={'a1':1 if binary else 0,'a2':0,'a3':0,'a4':0 if binary else 2,'a6':1 if binary else 2}
    g=(0,1) if binary else (5,1)
    c=Curve({'coefficients':{k:f.hex(v) for k,v in values.items()},'subgroup':{'order':'2' if binary else '19',
             'generator':{'x':f.hex(g[0]),'y':f.hex(g[1])},'r_torsion_rank':1}},f)
    return f,c


class ReferenceTests(unittest.TestCase):
    def test_prime_known_curve(self):
        f,c=toy('prime');g=(5,1)
        self.assertEqual(c.double(g),(6,3));self.assertEqual(c.mul(3,g),(10,6))
        self.assertIsNone(c.mul(19,g));self.assertEqual(c.mul(-1,g),(5,16))
        self.assertEqual(c.encode(g),'0305');self.assertEqual(c.decode('040501'),g)

    def test_binary_torsion_and_scalar_domain(self):
        f,c=toy('binary');p=(1,0);t=(0,1)
        self.assertTrue(c.on_curve(p));self.assertEqual(c.double(p),t)
        self.assertIsNone(c.double(t));self.assertEqual(c.neg(p),(1,1))
        self.assertEqual(c.mul(2,p),t) # r=2: reducing the scalar by r would be wrong.
        self.assertFalse(c.public_key_valid(p));self.assertTrue(c.public_key_valid(t))
        self.assertEqual(c.encode(p),'0201');self.assertEqual(c.encode((1,1)),'0301')
        self.assertEqual(c.decode('0200'),t)
        with self.assertRaises(ECError) as error:c.decode('0300')
        self.assertEqual(error.exception.code,'invalid_y_selector')

    def test_field_rejections_and_inverse(self):
        for kind in ('prime','binary'):
            f,c=toy(kind)
            for a in range(1,f.q):self.assertEqual(f.mul(a,f.inverse(a)),1)
            for invalid in (-1,f.q,True,1.5,'1'):
                for method in (f.inverse,f.square,f.sqrt):
                    with self.assertRaises(ECError):method(invalid)
                with self.assertRaises(ECError):f.mul(1,invalid)
            with self.assertRaises(ECError):f.inverse(0)

    def test_strict_encoding(self):
        f,c=toy('prime')
        for bad in ('','0000','0200ff','0511','0211','041101','040100'):
            with self.assertRaises(ECError):c.decode(bad)
        self.assertIsNone(c.decode('00'))
        zero=Curve({'coefficients':{k:f.hex(v) for k,v in {'a1':0,'a2':0,'a3':0,'a4':16,'a6':0}.items()},
          'subgroup':{'order':'2','generator':{'x':f.hex(0),'y':f.hex(0)},'r_torsion_rank':2}},f)
        with self.assertRaises(ECError) as error:zero.decode('0300')
        self.assertEqual(error.exception.code,'invalid_y_selector')

    def test_rank_two_membership(self):
        d=load_dataset();fields={f['id']:Field(f) for f in d['fields']}
        for cid in ('p61-E4-0','p64-E3-2','p64-E3-3'):
            r=next(c for c in d['curves'] if c['id']==cid);c=Curve(r,fields[r['field_id']])
            h=c.parse_point(r['named_points']['H']['point'])
            self.assertTrue(c.r_torsion(h));self.assertFalse(c.subgroup_contains(list(h)))
            self.assertTrue(c.subgroup_contains(list(c.G)));self.assertFalse(c.public_key_valid(None))

    def test_certificate_tamper(self):
        certs=copy.deepcopy(load_dataset()['certificates'])
        next(c for c in certs if c['n']=='3')['witnesses']['2']='1'
        with self.assertRaises((ValueError,AssertionError)):certify(certs)

    def test_full_release(self):
        result=verify(load_dataset())
        self.assertEqual(result['status'],'pass');self.assertEqual(result['curves'],56)


if __name__=='__main__':unittest.main()
