#!/usr/bin/env python3
"""Summarize recorded runs only; never manufacture certificate counts."""
from pathlib import Path
from collections import Counter,defaultdict
import json,statistics,sys
P=Path(__file__).resolve().parent.parent;L=P/'logs'
def read(s):
 p=L/(s+'.jsonl');return [json.loads(x) for x in p.open()] if p.exists() else []
def report(s,target):
 rows=read(s);ids={r['id']:r for r in rows};ct=Counter(r['status'] for r in ids.values());pr=[r for r in ids.values() if r['status']=='PROVED'];return dict(target=target,checked=len(ids),direct=sum(r.get('parts')==1 for r in pr),split=sum(r.get('parts',0)>1 for r in pr),unproved=ct['UNPROVED'],errors=ct['ERROR'],not_run=max(0,target-len(ids)),median_wall=statistics.median(r['seconds'] for r in ids.values()) if ids else None,counts=dict(ct))
A={'near_inner_full':582,'near_largest':24,'near_nonconvex':16,'near_random':1000,'far_sample':8,'far_E3_sample':3,'q_lg_2000':2000,'stdlib_near_candidate':1,'stdlib_candidate':1}
rows={s:report(s,n) for s,n in A.items()};N=set();F=set();Q=set()
for s in ('near_inner_full','near_largest','near_nonconvex','near_sample','stdlib_near_candidate'):
 N.update(r['id'] for r in read(s) if r['status']=='PROVED')
for s in ('far_sample','far_E3_sample'):F.update(r['id'] for r in read(s) if r['status']=='PROVED')
for s in ('q_lg_2000','q_sample','stdlib_candidate'):Q.update(r['id'] for r in read(s) if r['status']=='PROVED')
B={};
for name in ['B_easy','B_spot_fixed']:
 c=defaultdict(Counter)
 for x in read(name):
  for lab,ct in x['counts'].items():c[x['chart']+'/'+lab].update(ct)
 B[name]=dict(c)
add=Counter((x['status'],x.get('split')) for x in read('B_addendum_keys'));qf=Counter(x['status'] for x in read('q_fibre_sample'));anal=read('q_analytic');anat=Counter(x['method'] for x in anal);anids={x['id'] for x in anal};fibre_new={x['id'] for x in read('q_fibre_sample') if x['status']=='FIBRE_PROVED'}-anids
summary={'A':rows,'unique_proved':{'near':len(N),'far':len(F),'q_lg':len(Q)},'B':B,'addendum':{str(k):v for k,v in add.items()},'q_fibre':dict(qf),'q_analytic':dict(anat),'q_analytic_total':len(anids),'q_extra_fibre_proved':len(fibre_new)}
(P/'SUMMARY.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n');print(json.dumps(summary,indent=2,ensure_ascii=False))
