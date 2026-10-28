import numpy as np, sys
from scipy.spatial import cKDTree
from scipy.special import sph_harm_y
f=sys.argv[1]
lines=open(f).read().split('\n'); bx=[list(map(float,lines[i].split()[:2])) for i in (5,6,7)]
d=np.loadtxt(f,skiprows=9); r=d[:,1:4]; L=np.array([b[1]-b[0] for b in bx]); r=(r-np.array([b[0] for b in bx]))%L
t=cKDTree(r,boxsize=L); pairs=t.query_pairs(3.4,output_type='ndarray')
n=len(r); q=np.zeros((n,13),complex); cnt=np.zeros(n)
v=r[pairs[:,1]]-r[pairs[:,0]]; v-=L*np.round(v/L); rr=np.linalg.norm(v,axis=1)
th=np.arccos(v[:,2]/rr); ph=np.arctan2(v[:,1],v[:,0])
for m in range(-6,7):
    y=sph_harm_y(6,m,th,ph); yo=sph_harm_y(6,m,np.pi-th,ph+np.pi)
    np.add.at(q[:,m+6],pairs[:,0],y); np.add.at(q[:,m+6],pairs[:,1],yo)
np.add.at(cnt,pairs[:,0],1); np.add.at(cnt,pairs[:,1],1)
q/=np.maximum(cnt,1)[:,None]; q6=np.sqrt(4*np.pi/13*np.sum(abs(q)**2,1))
z=r[:,2]; bins=np.linspace(0,L[2],13); k=np.digitize(z,bins)
prof=[q6[k==i].mean() for i in range(1,13)]
print(f.split('/')[-1],'Q6 by z-slab:',np.round(prof,2),' frac solid-like (Q6>0.35):',round((q6>0.35).mean(),2))
qn=q/np.linalg.norm(q,axis=1)[:,None]
dot=np.real(np.sum(qn[pairs[:,0]]*np.conj(qn[pairs[:,1]]),1))
con=np.zeros(n); np.add.at(con,pairs[:,0],dot>0.5); np.add.at(con,pairs[:,1],dot>0.5)
sol=con>=7
print('  solid fraction by z-slab:',np.round([sol[k==i].mean() for i in range(1,13)],2),' total',round(sol.mean(),2))
