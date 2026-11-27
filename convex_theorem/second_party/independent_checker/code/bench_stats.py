#!/usr/bin/env python3
"""Summarize completed benchmark records only. Censored rows never counted as failures or successes."""
import argparse,collections,glob,json,math,pathlib,statistics
from driver import STAGES,canonical
TOTAL={'near':38280,'q_lg':66138,'q_other':118962,'far_std':293846,'far_E2':47804,'far_E3':35865,'B_std':934801,'B_E2':185823,'B_E3':226305}
def pct(a):return sorted(a)[max(0,math.ceil(.95*len(a))-1)] if a else None
def med(a):return statistics.median(a) if a else None
def mean(a):return statistics.mean(a) if a else None
def main():
 p=argparse.ArgumentParser();p.add_argument('bench',nargs='+',type=pathlib.Path);p.add_argument('--out',type=pathlib.Path,required=True);a=p.parse_args();groups=collections.defaultdict(dict)
 for d in a.bench:
  for f in d.glob('*.jsonl'):
   for l in f.open():
    r=json.loads(l);groups[f.stem][r['id']]=r
 out={}
 for g,rr in sorted(groups.items()):
  rs=list(rr.values());ct=collections.Counter(r['status'] for r in rs);t=[r['seconds'] for r in rs];w=[len(canonical(r['witness']).encode()) for r in rs if 'witness'in r];rows=[len(canonical({k:v for k,v in r.items() if k not in ('task','old_status')}).encode())+1 for r in rs];o={'n':len(rs),'counts':dict(ct),'median_s':med(t),'p95_s':pct(t),'mean_s':mean(t),'max_s':max(t),'witness_median_bytes':med(w),'witness_p95_bytes':pct(w),'mean_row_bytes':mean(rows),'peak_worker_RSS_MiB':max(r['rss_kib'] for r in rs)/1024}
  if g.startswith('bench_'):
   stage=g[6:];n=TOTAL[stage];jobs=10 if stage in ('far_E2','far_E3') else 20;o.update(target=200,full=n,configured_jobs=jobs,ideal20_hours=mean(t)*n/20/3600,configured_hours=mean(t)*n/jobs/3600,raw_GiB=mean(rows)*n/2**30)
  elif g.startswith('focus_'):
   old=[r for r in rs if r.get('old_status')=='UNPROVED'];rc=collections.Counter();sp=0
   for r in old:
    if 'witness' not in r:continue
    def rec(w):
     if 'leaf' in w:rc[w['leaf']['route']]+=1
     else:rec(w['left']);rec(w['right'])
    rec(r['witness']);sp+=r['parts']>1
   o.update(old_fail_attempted=len(old),old_fail_new_counts=dict(collections.Counter(r['status'] for r in old)),old_fail_split=sp,old_fail_terminal_routes=dict(rc))
  out[g]=o
 a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':main()
