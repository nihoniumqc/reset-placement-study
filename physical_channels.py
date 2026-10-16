"""Physical Pauli-channel signatures by backward circuit propagation."""
import itertools

def physical_channels(c):
 c=c.flattened(); rec=[0]*c.num_measurements; starts={};m=d=0
 for off,ins in enumerate(c):
  if ins.name in ['M','MX','MR','MRX']:
   starts[off]=m;m+=len(ins.targets_copy())
  elif ins.name in ['DETECTOR','OBSERVABLE_INCLUDE']:
   bit=d if ins.name=='DETECTOR' else c.num_detectors+int(ins.gate_args_copy()[0])
   for t in ins.targets_copy():
    assert t.is_measurement_record_target
    rec[m+t.value]^=1<<bit
   if ins.name=='DETECTOR':d+=1
 assert m==c.num_measurements and d==c.num_detectors
 x=[0]*c.num_qubits;z=x.copy();out=[]
 for off in reversed(range(len(c))):
  ins=c[off];name=ins.name;ts=ins.targets_copy();args=ins.gate_args_copy()
  if name in ['DETECTOR','OBSERVABLE_INCLUDE','QUBIT_COORDS','SHIFT_COORDS','TICK']:continue
  assert all(t.is_qubit_target for t in ts),(name,ts)
  qs=[t.value for t in ts]
  if name in ['DEPOLARIZE1','DEPOLARIZE2','X_ERROR','Y_ERROR','Z_ERROR','PAULI_CHANNEL_1','PAULI_CHANNEL_2']:
   if not (sum(args) if name.startswith('PAULI_CHANNEL') else args[0]):continue
   stride=2 if name in ['DEPOLARIZE2','PAULI_CHANNEL_2'] else 1
   strings=list(itertools.product('IXYZ',repeat=stride))[1:] if name.startswith(('DEPOLARIZE','PAULI_CHANNEL')) else [(name[0],)]
   for start in range(0,len(qs),stride):
    ev=[]
    for event_index,letters in enumerate(strings):
     v=0
     for q,p in zip(qs[start:start+stride],letters):
      if p in 'XY':v^=x[q]
      if p in 'ZY':v^=z[q]
     ev.append((v,args[event_index] if name.startswith('PAULI_CHANNEL') else args[0]/len(strings)))
    out.append({'offset':off,'start':start,'gate':name,'p':sum(args) if name.startswith('PAULI_CHANNEL') else args[0],'events':ev})
  elif name in ['M','MX','MR','MRX']:
   for j,q in enumerate(qs):
    if args and args[0]:out.append({'offset':off,'start':j,'gate':name,'p':args[0],'events':[(rec[starts[off]+j],args[0])]})
    if name in ['MR','MRX']:x[q]=z[q]=0
    if name in ['M','MR']:x[q]^=rec[starts[off]+j]
    else:z[q]^=rec[starts[off]+j]
  elif name in ['R','RX']:
   for q in qs:x[q]=z[q]=0
  elif name=='H':
   for q in qs:x[q],z[q]=z[q],x[q]
  elif name in ['S','S_DAG']:
   for q in qs:x[q]^=z[q]
  elif name=='CX':
   for j in reversed(range(0,len(qs),2)):
    a,b=qs[j:j+2];x[a]^=x[b];z[b]^=z[a]
  elif name not in ['I','X','Y','Z']:raise ValueError(name)
 return sorted(out,key=lambda r:(r['offset'],r['start']))
