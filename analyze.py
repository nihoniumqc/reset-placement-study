"""Reconstruct paired effects and simultaneous bounds from selected aggregate data."""
from pathlib import Path
import json
import hashlib
import stim
import numpy as np
from scipy.stats import beta
ROOT=Path(__file__).resolve().parent

def interval(k,n,family):
 alpha=.05/family
 return (0. if k==0 else float(beta.ppf(alpha/2,k,n-k+1)),1. if k==n else float(beta.ppf(1-alpha/2,k+1,n-k)))

def check_dense_results():
 data=json.loads((ROOT/'data/dense_results.json').read_text())
 rows=data['rows'];family=data['family_marginals']
 assert family==len(rows)*10==270
 conditions={c['id']:c for c in json.loads((ROOT/'data/conditions.json').read_text())}
 for r in rows:
  c=conditions[r['condition']]
  assert r['circuit_sha256']==c['circuit_sha256']
  n=r['shots'];a,b,repair,harm=r['counts'];tc=np.array(r['target_counts'])
  assert tc.shape==(3,2) and a-b==repair-harm
  assert repair==tc[:,1].sum() and harm==tc[:,0].sum()
  marginal=np.array([interval(k,n,family) for k in r['counts']+tc.ravel().tolist()])
  ci=[marginal[2,0]-marginal[3,1],marginal[2,1]-marginal[3,0]]
  assert np.allclose(marginal[:4],r['display_marginals'],atol=1e-14,rtol=0)
  assert np.allclose(ci,r['display_interval'],atol=1e-14,rtol=0)
  assert abs((repair-harm)/n-r['effect'])<1e-14
  for j,t in enumerate(r['syndromes']):
   zero,one=marginal[4+2*j:6+2*j]
   assert np.allclose([one[0]-zero[1],one[1]-zero[0]],t['interval'],atol=1e-14,rtol=0)
   assert abs((tc[j,1]-tc[j,0])/n-t['effect'])<1e-14
 assert len({r['scale'] for r in rows if r['kind']=='noise'})==17
 assert len({r['factor'] for r in rows if r['kind']=='attenuation'})==9
 for c in conditions.values():
  assert hashlib.sha256(str(stim.Circuit((ROOT/c['circuit']).read_text())).encode()).hexdigest()==c['circuit_sha256']
 target=next(r for r in rows if r['kind']=='attenuation' and r['factor']==0)
 control=next(r for r in rows if r['kind']=='matched')
 baseline=next(r for r in rows if r['kind']=='attenuation' and r['factor']==1)
 for label,left,right in [('target - matched',target,control),('target - baseline',target,baseline),('matched - baseline',control,baseline)]:
  value=left['effect']-right['effect']
  ci=[left['display_interval'][0]-right['display_interval'][1],left['display_interval'][1]-right['display_interval'][0]]
  print(f'{label:24s} {value*1e6:9.3f} ppm [{ci[0]*1e6:9.3f}, {ci[1]*1e6:9.3f}]')
 print(f"{len(rows)} dense-sweep effects and all {family} displayed marginal intervals reproduced; {len(conditions)} circuit hashes verified.")

def main():
 rows=json.loads((ROOT/'data/paired_results.json').read_text())
 for r in rows:
  a,b,repair,harm=r['counts'];n=r['shots'];assert a-b==repair-harm
  estimate=(repair-harm)/n;lo1,hi1=interval(repair,n,r['family_marginals']);lo0,hi0=interval(harm,n,r['family_marginals']);ci=[lo1-hi0,hi1-lo0]
  assert np.allclose(ci,r['interval'],atol=1e-14,rtol=0) and abs(estimate-r['effect'])<1e-14
  print(f"{r['condition']:16s} {r['policy']:10s} {estimate*1e6:9.3f} ppm [{ci[0]*1e6:9.3f}, {ci[1]*1e6:9.3f}]")
 cor=json.loads((ROOT/'data/correlated_results.json').read_text())
 for r in cor:
  for v in r['comparisons']:
   assert r['counts'][v['left']]-r['counts'][v['right']]==v['repairs']-v['harms']
   lo1,hi1=interval(v['repairs'],r['shots'],r['family_marginals']);lo0,hi0=interval(v['harms'],r['shots'],r['family_marginals'])
   assert np.allclose([lo1-hi0,hi1-lo0],v['interval'],atol=1e-14,rtol=0)
 check_dense_results()
 print(f'{len(rows)} paired effects and {sum(len(r["comparisons"]) for r in cor)} correlated contrasts reproduced.')
if __name__=='__main__':main()
