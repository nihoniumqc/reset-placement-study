"""All local fault orders convolved with at most two background faults.

This is an unnormalized approximation, not a certified near-ML decoder.
"""
from pathlib import Path
from collections import defaultdict
import json,math,argparse
import numpy as np
import stim,pymatching
from physical_channels import physical_channels
ROOT=Path(__file__).resolve().parent
SUPPORT=[60,63,64,66,69,70,74,77,79,80,83,89,93,94]
EXTENDED=[76,79,80,82,85,86,90,93,95,96,99,105,109,110]

def calculate(condition):
 c=stim.Circuit((ROOT/condition['circuit']).read_text());nd=c.num_detectors
 support=EXTENDED if condition['experiment']=='boundary' else SUPPORT
 outside=((1<<nd)-1)^sum(1<<d for d in support);dim=1<<(len(support)+1);lbit=dim//2;idx=np.arange(dim)
 def project(v):return sum(((v>>d)&1)<<j for j,d in enumerate(support))|((v>>nd)<<len(support))
 chs=physical_channels(c);lc=[ch for ch in chs if all(v&outside==0 for v,p in ch['events'])];bg=[ch for ch in chs if any(v&outside for v,p in ch['events'])]
 local=np.zeros(dim);local[0]=1
 for ch in lc:
  probs=defaultdict(float);probs[0]=1-ch['p']
  for v,p in ch['events']:probs[project(v)]+=p
  local=sum(p*local[idx^v] for v,p in probs.items())
 assert abs(local.sum()-1)<1e-12
 p0=math.prod(1-ch['p'] for ch in bg);groups=defaultdict(lambda:defaultdict(float));same=np.zeros(dim);q0,q1,q2=1.,0.,0.
 for ch in bg:
  p=ch['p'];q0,q1,q2=q0*(1-p),q1*(1-p)+q0*p,q2*(1-p)+q1*p
  ev=[(v&outside,project(v),p/(1-ch['p'])) for v,p in ch['events']]
  for o,i,w in ev:groups[o][i]+=w
  for o,i,w in ev:
   for u,j,t in ev:
    if o==u:same[i^j]+=w*t
 b0=np.zeros(dim);b0[0]=p0;b1=np.zeros(dim);b2=-same
 for v,w in groups[0].items():b1[v]=p0*w
 for g in groups.values():
  for i,w in g.items():
   for j,t in g.items():b2[i^j]+=w*t
 b2*=p0/2;assert b2.min()>-1e-13;b2=np.maximum(b2,0)
 total=b0+b1+b2;codes=[project(sum(1<<d for d in t)) for t in condition['targets']]
 qs=[[float(local[idx^t^(label*lbit)]@total) for label in [0,1]] for t in codes]
 labels=[int(q[1]>q[0]) for q in qs]
 m=pymatching.Matching.from_detector_error_model(c.detector_error_model(decompose_errors=True,approximate_disjoint_errors=True))
 ds=np.zeros((len(codes),nd),dtype=np.uint8)
 for i,t in enumerate(condition['targets']):ds[i,t]=1
 baseline=m.decode_batch(ds)[:,0].tolist()
 return {'condition':condition['id'],'q_approx':qs,'predicted_labels':labels,'ordinary_labels':baseline,'omitted_background_mass':1-q0-q1-q2,'label_one_effect':sum(q[1]-q[b] for q,b in zip(qs,baseline)),'predicted_effect':sum(q[l]-q[b] for q,l,b in zip(qs,labels,baseline))}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('condition');a=ap.parse_args()
 configs={r['id']:r for r in json.loads((ROOT/'data/conditions.json').read_text())}
 print(json.dumps(calculate(configs[a.condition]),indent=2))
