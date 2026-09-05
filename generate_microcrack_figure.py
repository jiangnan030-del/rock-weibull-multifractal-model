"""Generate an uploaded SVG overview of Monte Carlo crack networks."""
from pathlib import Path
import numpy as np
from microcrack_multifractal import generate_microcrack_network
SPACINGS=(.9,.7,.5,.3,.1); CUBE=10.
def project(p):
 p=np.asarray(p); return np.column_stack(((p[:,0]-p[:,1])*.866,(p[:,0]+p[:,1])*.5-p[:,2]))
def panel(spacing,ox,oy,size,seed):
 n=generate_microcrack_network(CUBE,spacing,24,.18,seed); corners=np.array([[x,y,z] for x in (0,CUBE) for y in (0,CUBE) for z in (0,CUBE)],float); pc=project(corners); lo=pc.min(0); span=np.ptp(pc,axis=0); scale=size/max(span)
 def xy(p):
  q=(project(np.asarray(p).reshape(-1,3))-lo)*scale; q[:,0]+=ox+(size-span[0]*scale)/2; q[:,1]+=oy+(size-span[1]*scale)/2; return q
 out=[]
 for i,a in enumerate(corners):
  for j,b in enumerate(corners):
   if j>i and np.count_nonzero(a!=b)==1:
    q=xy([a,b]); out.append(f'<line x1="{q[0,0]:.1f}" y1="{q[0,1]:.1f}" x2="{q[1,0]:.1f}" y2="{q[1,1]:.1f}" class="box"/>')
 for cid in range(len(n.centers)):
  p=n.points[n.crack_ids==cid]
  if len(p)<2: continue
  c=p-p.mean(0); *_,vh=np.linalg.svd(c,full_matrices=False); axis=vh[0]; t=c@axis; q=xy([p.mean(0)+axis*t.min(),p.mean(0)+axis*t.max()]); out.append(f'<line x1="{q[0,0]:.1f}" y1="{q[0,1]:.1f}" x2="{q[1,0]:.1f}" y2="{q[1,1]:.1f}" class="crack"/>')
 out.append(f'<text x="{ox+size/2:.1f}" y="{oy+size+25:.1f}" text-anchor="middle">spacing {spacing:g} cm</text>'); return '\n'.join(out)
def main():
 positions=((45,30),(410,30),(775,30),(225,390),(590,390)); body='\n'.join(panel(s,*pos,280,2026+i) for i,(s,pos) in enumerate(zip(SPACINGS,positions))); svg=f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1100 720"><style>line{{vector-effect:non-scaling-stroke}}.box{{stroke:#555;stroke-width:1}}.crack{{stroke:#111;stroke-width:.8;opacity:.82}}text{{font:18px Arial}}</style><rect width="100%" height="100%" fill="white"/>{body}</svg>'''; p=Path('figures/figure1_microcrack_networks.svg'); p.parent.mkdir(exist_ok=True); p.write_text(svg)
if __name__=='__main__': main()
