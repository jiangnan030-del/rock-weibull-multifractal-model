"""3D Monte Carlo crack network and Chhabra-Jensen multifractal spectrum."""
from __future__ import annotations
import argparse,csv,json,math
from dataclasses import asdict,dataclass
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
@dataclass(frozen=True)
class CrackSet: set_id:int; direction:float; direction_var:float; dip:float; dip_var:float; length:float; length_var:float
TABLE1_CRACK_SETS=(CrackSet(1,2.7,145.9,75.7,359.9,4.,37.),CrackSet(2,340.7,419.3,34.,106.8,6.29,8.75))
@dataclass
class CrackNetwork: points:np.ndarray; crack_ids:np.ndarray; set_ids:np.ndarray; centers:np.ndarray; cube_size_cm:float; target_spacing_cm:float; seed:int
@dataclass(frozen=True)
class ScaleStatistics: division:int; r:float; occupied_boxes:int; probabilities:np.ndarray
@dataclass(frozen=True)
class LinearFit: slope:float; intercept:float; correlation:float; r_squared:float
@dataclass(frozen=True)
class SpectrumPoint:
 q:float; alpha:float; f:float; start_index:int; end_index:int; min_division:int; max_division:int; correlation_alpha:float; correlation_f:float; r2_alpha:float; r2_f:float; threshold_passed:bool
def _basis(direction,dip):
 a,d=np.radians([direction,dip]); n=np.array([np.sin(d)*np.sin(a),np.sin(d)*np.cos(a),np.cos(d)]); h=np.array([0.,0.,1.]) if abs(n@np.array([0.,0.,1.]))<=.95 else np.array([1.,0.,0.]); u=np.cross(n,h); u/=np.linalg.norm(u); return u,np.cross(n,u)
def generate_microcrack_network(cube_size_cm=10.,target_spacing_cm=.5,points_per_crack=300,cracks_per_spacing_per_set=1.,seed=2026,max_trace_fraction=.9):
 if cube_size_cm<=0 or target_spacing_cm<=0: raise ValueError('cube size and spacing must be positive')
 if points_per_crack<8: raise ValueError('points_per_crack must be >= 8')
 rng=np.random.default_rng(seed); count=max(1,round(cube_size_cm/target_spacing_cm*cracks_per_spacing_per_set)); pts=[]; cids=[]; sids=[]; centers=[]; cid=0
 for s in TABLE1_CRACK_SETS:
  directions=rng.normal(s.direction,np.sqrt(s.direction_var),count)%360; dips=np.abs((rng.normal(s.dip,np.sqrt(s.dip_var),count)+90)%180-90); lengths=np.clip(rng.exponential(s.length,count),cube_size_cm*.01,cube_size_cm*max_trace_fraction); cs=rng.uniform(0,cube_size_cm,(count,3))
  for direction,dip,length,center in zip(directions,dips,lengths,cs,strict=True):
   u,v=_basis(direction,dip); theta=rng.uniform(0,2*np.pi,points_per_crack); radius=.5*length*np.sqrt(rng.random(points_per_crack)); p=center+radius[:,None]*(np.cos(theta)[:,None]*u+np.sin(theta)[:,None]*v); p=p[np.all((p>=0)&(p<=cube_size_cm),axis=1)]; p=p if len(p) else center[None,:]; pts.append(p); cids.append(np.full(len(p),cid)); sids.append(np.full(len(p),s.set_id)); centers.append(center); cid+=1
 return CrackNetwork(np.vstack(pts),np.concatenate(cids),np.concatenate(sids),np.vstack(centers),cube_size_cm,target_spacing_cm,seed)
def box_count_statistics(points,cube_size_cm,divisions):
 if np.asarray(points).ndim!=2 or np.asarray(points).shape[1]!=3 or len(points)==0: raise ValueError('points must be nonempty (n,3)')
 p=np.clip(points/cube_size_cm,0,np.nextafter(1.,0.)); out=[]
 for n in sorted(set(divisions)):
  if n<2: raise ValueError('box divisions must be >= 2')
  idx=np.floor(p*n).astype(int); linear=idx[:,0]*n*n+idx[:,1]*n+idx[:,2]; _,counts=np.unique(linear,return_counts=True); prob=counts/counts.sum(); out.append(ScaleStatistics(n,1/n,len(counts),prob))
 if len(out)<3: raise ValueError('at least 3 scales required')
 return out
def _fit(x,y):
 if np.allclose(y,y[0]): return LinearFit(0.,float(y[0]),1.,1.)
 slope,intercept=np.polyfit(x,y,1); corr=float(np.corrcoef(x,y)[0,1]); return LinearFit(float(slope),float(intercept),corr,corr*corr)
def chhabra_jensen_scale_values(scales,q):
 xs=[]; ya=[]; yf=[]
 for s in scales:
  lp=np.log(s.probabilities); lw=q*lp; w=np.exp(lw-lw.max()); mu=w/w.sum(); xs.append(np.log(s.r)); ya.append(np.sum(mu*lp)); yf.append(np.sum(mu*np.log(mu)))
 return np.array(xs),np.array(ya),np.array(yf)
def select_scaling_window(x,y_alpha,y_f,correlation_threshold=.95,min_points=3):
 if not(len(x)==len(y_alpha)==len(y_f)): raise ValueError('arrays must have equal length')
 if min_points<3 or min_points>len(x): raise ValueError('invalid min_points')
 passed=[]; physical=[]; all_windows=[]
 for start in range(len(x)-min_points+1):
  for end in range(start+min_points,len(x)+1):
   a,f=_fit(x[start:end],y_alpha[start:end]),_fit(x[start:end],y_f[start:end]); score=(end-start,min(a.r_squared,f.r_squared),(a.r_squared+f.r_squared)/2,-start,start,end,a,f); all_windows.append(score); valid=a.slope>=0 and 0<=f.slope<=3
   if valid: physical.append(score)
   if valid and abs(a.correlation)>=correlation_threshold and abs(f.correlation)>=correlation_threshold: passed.append(score)
 candidates=passed or physical or all_windows
 if passed: chosen=max(candidates,key=lambda z:z[:4]); ok=True
 else: chosen=max(candidates,key=lambda z:(z[1],z[2],z[0],z[3])); ok=False
 *_,start,end,a,f=chosen; return start,end,a,f,ok
def compute_multifractal_spectrum(scales,q_values,correlation_threshold=.95,min_points=3):
 out=[]
 for q in q_values:
  x,ya,yf=chhabra_jensen_scale_values(scales,float(q)); start,end,a,f,ok=select_scaling_window(x,ya,yf,correlation_threshold,min_points); ds=[s.division for s in scales[start:end]]; out.append(SpectrumPoint(float(q),a.slope,f.slope,start,end,min(ds),max(ds),a.correlation,f.correlation,a.r_squared,f.r_squared,ok))
 return out
def _csv(path,header,rows):
 path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('w',newline='',encoding='utf-8') as h: w=csv.writer(h); w.writerow(header); w.writerows(rows)
def save_outputs(network,scales,spectrum,out):
 out.mkdir(parents=True,exist_ok=True); _csv(out/'microcrack_points.csv',['x_cm','y_cm','z_cm','crack_id','crack_set'],((*p,int(c),int(s)) for p,c,s in zip(network.points,network.crack_ids,network.set_ids,strict=True))); _csv(out/'box_count_scales.csv',['division','relative_box_size_r','occupied_boxes'],((s.division,s.r,s.occupied_boxes) for s in scales))
 with (out/'multifractal_spectrum.csv').open('w',newline='',encoding='utf-8') as h: w=csv.DictWriter(h,fieldnames=asdict(spectrum[0]).keys()); w.writeheader(); w.writerows(asdict(s) for s in spectrum)
 rng=np.random.default_rng(network.seed+1); ids=np.arange(len(network.points)); ids=rng.choice(ids,20000,replace=False) if len(ids)>20000 else ids; fig=plt.figure(figsize=(7,6)); ax=fig.add_subplot(111,projection='3d')
 for sid,color in ((1,'tab:blue'),(2,'tab:orange')):
  p=network.points[ids][network.set_ids[ids]==sid]; ax.scatter(*p.T,s=1,alpha=.25,color=color,label=f'set {sid}')
 ax.set(xlabel='x (cm)',ylabel='y (cm)',zlabel='z (cm)'); ax.legend(); fig.tight_layout(); fig.savefig(out/'microcrack_network_3d.png',dpi=180); plt.close(fig)
 alpha=np.array([s.alpha for s in spectrum]); fv=np.array([s.f for s in spectrum]); ok=np.array([s.threshold_passed for s in spectrum]); fig,ax=plt.subplots(figsize=(6,5)); ax.plot(alpha,fv,'-',color='.6'); ax.scatter(alpha[ok],fv[ok],label='passed')
 if np.any(~ok): ax.scatter(alpha[~ok],fv[~ok],marker='x',color='red',label='fallback')
 ax.set(xlabel=r'$\alpha(q)$',ylabel=r'$f(q)$'); ax.grid(alpha=.3); ax.legend(); fig.tight_layout(); fig.savefig(out/'multifractal_spectrum.png',dpi=180); plt.close(fig)
def main():
 p=argparse.ArgumentParser(); p.add_argument('--cube-size',type=float,default=10.); p.add_argument('--spacing',type=float,default=.5); p.add_argument('--points-per-crack',type=int,default=300); p.add_argument('--crack-density-factor',type=float,default=1.); p.add_argument('--seed',type=int,default=2026); p.add_argument('--divisions',default='4,5,6,8,10,12,16,20'); p.add_argument('--q-min',type=float,default=-5); p.add_argument('--q-max',type=float,default=5); p.add_argument('--q-step',type=float,default=.5); p.add_argument('--correlation-threshold',type=float,default=.95); p.add_argument('--min-window-points',type=int,default=4); p.add_argument('--outdir'); a=p.parse_args(); divisions=[int(x) for x in a.divisions.split(',')]; q=np.arange(a.q_min,a.q_max+a.q_step/2,a.q_step); out=Path(a.outdir) if a.outdir else Path('outputs')/f'multifractal_spacing_{a.spacing:g}'; network=generate_microcrack_network(a.cube_size,a.spacing,a.points_per_crack,a.crack_density_factor,a.seed); scales=box_count_statistics(network.points,network.cube_size_cm,divisions); spectrum=compute_multifractal_spectrum(scales,q,a.correlation_threshold,a.min_window_points); save_outputs(network,scales,spectrum,out); q0=min(spectrum,key=lambda s:abs(s.q)); summary={'cube_size_cm':a.cube_size,'target_spacing_cm':a.spacing,'seed':a.seed,'crack_count':len(network.centers),'sampled_point_count':len(network.points),'box_divisions':divisions,'capacity_dimension_D0':q0.f,'relative_dimension_lambda':q0.f/3,'spectrum_width':max(s.alpha for s in spectrum)-min(s.alpha for s in spectrum),'all_q_windows_passed':all(s.threshold_passed for s in spectrum),'table1_crack_sets':[asdict(s) for s in TABLE1_CRACK_SETS]}; (out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8'); print(json.dumps(summary,ensure_ascii=False,indent=2))
if __name__=='__main__': main()
