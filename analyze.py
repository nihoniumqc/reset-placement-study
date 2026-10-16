"""Reconstruct paired effects and simultaneous bounds from selected aggregate data."""
from pathlib import Path
import json
import numpy as np
from scipy.stats import beta
ROOT=Path(__file__).resolve().parent

def interval(k,n,family):
 alpha=.05/family
 return (0. if k==0 else float(beta.ppf(alpha/2,k,n-k+1)),1. if k==n else float(beta.ppf(1-alpha/2,k+1,n-k)))

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
 print(f'{len(rows)} paired effects and {sum(len(r["comparisons"]) for r in cor)} correlated contrasts reproduced.')
if __name__=='__main__':main()
