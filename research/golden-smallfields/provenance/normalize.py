import csv
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
p=json.loads((ROOT/'prime_data.json').read_text())
b=json.loads((ROOT/'binary_data.json').read_text())
out={'fields':[], 'curves':[]}
for field in p['fields']:
    out['fields'].append({'id':field['id'],'kind':'prime','q':field['p'],'bits':field['bits'],'form':field['form'],'g':field['primitive_root']})
    for family in field['families']:
        for row in family['twists']:
            family_name=family['family']
            idx=row['index']
            c={'id':f"{field['id']}-{family_name}-{idx}", 'field':field['id'], 'family':family_name, 'kind':'prime',
               'q':field['p'],'t':row['trace'],'N':row['order'],'r':row['largest_prime_factor'],'h':row['cofactor_for_largest_prime'],
               'k':row['embedding_degree'],'factors':row['order_factors'],'group':row['group_invariant_factors'],
               'Dpi':row['Frobenius_discriminant'],'DK':row['quadratic_field_discriminant'],'fpi':row['Frobenius_order_conductor'],
               'Dend':row['endomorphism_discriminant_over_Fp'],'ordinary':family['status']=='ordinary',
               'A':row['A'],'B':row['B'],'j':row.get('j_mod_p',family['j_mod_p']), 'delta':row['curve_discriminant_mod_p'],
               'index':idx,'g':field['primitive_root'], 'twist_id':f"{field['id']}-{family_name}-{row['quadratic_twist_index']}",
               'endo_proof':row.get('endomorphism_proof','Known special-CM model or supersingular rational 2-torsion test.')}
            out['curves'].append(c)
s=p['searched_example']
for row in s['twists']:
    out['curves'].append({'id':f"p83-S-{row['index']}",'field':'p83','family':'S','kind':'prime','q':s['p'],'t':row['trace'],
      'N':row['order'],'r':row['largest_prime_factor'],'h':row['cofactor_for_largest_prime'],'k':row['embedding_degree'],
      'factors':row['order_factors'],'group':row['group_invariant_factors'],'Dpi':row['Frobenius_discriminant'],
      'DK':row['quadratic_field_discriminant'],'fpi':row['Frobenius_order_conductor'],'Dend':row['endomorphism_discriminant_over_Fp'],
      'ordinary':True,'A':row['A'],'B':row['B'],'j':row['j_mod_p'],'delta':row['curve_discriminant_mod_p'],'index':row['index'],
      'twist_id':f"p83-S-{1-row['index']}",'endo_proof':row['endomorphism_proof']})
for field in b['fields']:
    m=field['degree'];out['fields'].append({'id':f'b{m}','kind':'binary','q':2**m,'bits':m,'exponents':field['exponents']})
for key, family in [('curves','K'),('generic_curves','B'),('searched_examples','S')]:
    for row in b[key]:
        out['curves'].append({'id':row['label'],'field':f"b{row['field_degree']}",'family':family,'kind':'binary',
          'q':row['q'],'t':row['trace'],'N':row['group_order'],'r':row['largest_prime_subgroup'],'h':row['cofactor'],
          'k':row['embedding_degree'],'factors':row['group_order_factors'],'group':row['group_invariant_factors'],
          'Dpi':row['frobenius_discriminant'],'DK':row['CM_field_discriminant'],'fpi':row['frobenius_order_conductor'],
          'Dend':row['endomorphism_ring_discriminant'],'ordinary':True,'a':row['a'],'b':row['b'],
          'j':row['j_invariant'],'delta':row['curve_discriminant_in_field'],'twist_id':row['twist_label'],
          'frob_classno':row.get('frobenius_order_class_number'), 'b_exponent':row.get('b_exponent',0)})
(ROOT/'normalized_data.json').write_text(json.dumps(out,indent=2)+'\n')
columns=['id','field','family','q','t','N','r','h','k','Dpi','DK','fpi','Dend','ordinary']
with (ROOT/'curve_invariants.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,columns,extrasaction='ignore');w.writeheader();w.writerows(out['curves'])
print('Normalized',len(out['fields']),'fields and',len(out['curves']),'curve records.')
