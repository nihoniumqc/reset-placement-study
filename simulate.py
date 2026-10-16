"""Replay selected stored circuits with ordinary or correlated matching.

Example: python simulate.py mixed_00 --shots 10000
Use --original-shots only when intentionally replaying the full budget.
"""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import stim,pymatching
ROOT=Path(__file__).resolve().parent

def simulate(condition,shots,correlated=False):
 c=stim.Circuit((ROOT/condition['circuit']).read_text())
 assert hashlib.sha256(str(c).encode()).hexdigest()==condition['circuit_sha256']
 dem=c.detector_error_model(decompose_errors=True,approximate_disjoint_errors=True)
 m=pymatching.Matching.from_detector_error_model(dem,enable_correlations=correlated)
 counts={p:np.zeros(4,dtype=np.int64) for p in condition['policies']}
 for chunk,start in enumerate(range(0,shots,condition['chunk_size'])):
  n=min(condition['chunk_size'],shots-start)
  ds,obs=c.compile_detector_sampler(seed=condition['seed_base']+chunk).sample(n,separate_observables=True)
  pred=m.decode_batch(ds,enable_correlations=correlated)[:,0];obs=obs[:,0]
  hs=np.column_stack([(ds.sum(axis=1)==len(t)) & np.all(ds[:,t].astype(bool),axis=1) for t in condition['targets']]).astype(bool)
  old=pred!=obs
  for policy,labels in condition['policies'].items():
   alt=pred.copy()
   for j,label in enumerate(labels):alt[hs[:,j]]=label
   new=alt!=obs
   counts[policy]+=np.array([old.sum(),new.sum(),(old&~new).sum(),(~old&new).sum()])
 return {'condition':condition['id'],'shots':shots,'correlated':correlated,'policies':{p:{'counts':v.tolist(),'effect':float((v[2]-v[3])/shots)} for p,v in counts.items()}}

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('condition');ap.add_argument('--shots',type=int,default=10000);ap.add_argument('--original-shots',action='store_true');ap.add_argument('--correlated',action='store_true');a=ap.parse_args()
 cs={c['id']:c for c in json.loads((ROOT/'data/conditions.json').read_text())};c=cs[a.condition];n=c['original_shots'] if a.original_shots else a.shots
 if n<1:ap.error('shots must be positive')
 print(json.dumps(simulate(c,n,a.correlated),indent=2))
