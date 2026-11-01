#!/usr/bin/env python3
"""Independent MPS FHM verifier. It NEVER imports producer code.
Use mpmath interval arithmetic and a Bessel Taylor series with explicit tails.
Untrusted input is read as data only. Runs for selected 3 exported MPS traces.
"""
import json,math,sys,time
from mpmath import iv
iv.dps=35
P=iv.pi

def Q(f):
    x=float(f)
    a,b=x.as_integer_ratio()
    return iv.mpf(a)/iv.mpf(b)

def high(x):return math.nextafter(float(x.b),math.inf)
def low(x):return math.nextafter(float(x.a),-math.inf)
def abmax(x):return max(math.nextafter(abs(float(x.a)),math.inf),math.nextafter(abs(float(x.b)),math.inf))
def bu(x):return iv.mpf([str(x.a)[1:].split(',')[0],str(x.b)[1:].split(',')[0]])
def dot(x,y):return x[0]*y[0]+x[1]*y[1]
def cross(x,y):return x[0]*y[1]-x[1]*y[0]
def diff(a,b):return [a[0]-b[0],a[1]-b[1]]
def add(a,b):return [a[0]+b[0],a[1]+b[1]]
def scal(c,a):return [c*a[0],c*a[1]]

def bessel_series(nu,z,gamma_inv,N=38):
    """Return enclosing (J_nu, J'_nu, J''_nu).
    For n>=N, term ratios bounded by q, including derivative weights explicitly.
    """
    zz=z/2
    term=iv.exp(nu*iv.ln(zz))*gamma_inv
    summ=iv.mpf('0');ds=iv.mpf('0');dds=iv.mpf('0')
    for j in range(N):
        t=nu+2*j
        summ+=term;ds+=term*t/z;dds+=term*t*(t-1)/(z*z)
        term= -term*z*z/(4*(j+1)*(nu+j+1))
    q=(abmax(z)**2)/(4*(N+1)*(low(nu)+N+1))
    if q>=1:raise ValueError('Bessel tail ratio >=1')
    # For t=N term and following, series of |term| q^l and factors up to t0+2l.
    next_abs=abmax(term); zlow=low(z); t0=max(0,high(nu))+2*N
    if zlow<=0:raise ValueError('Bessel argument contains zero')
    geometric=1/(1-q)
    m1=q/(1-q)**2
    m2=q*(1+q)/(1-q)**3
    err0=next_abs*geometric
    err1=next_abs/zlow*(t0*geometric+2*m1)
    err2=next_abs/zlow**2*((t0*t0+t0)*geometric+ (4*t0+2)*m1+4*m2)
    # outward binary64 error margins on numerical tail radius
    err0=math.nextafter(err0,math.inf)*1.00000000001
    err1=math.nextafter(err1,math.inf)*1.00000000001
    err2=math.nextafter(err2,math.inf)*1.00000000001
    return summ+iv.mpf([-err0,err0]),ds+iv.mpf([-err1,err1]),dds+iv.mpf([-err2,err2])

def setup(data):
    V=[[Q(a) for a in row] for row in data['vertices']]
    frames=[]; coeffs=[Q(x) for x in data['coefficients']]
    k=iv.sqrt(Q(data['lam_tilde']))
    for j in range(5):
        e=diff(V[(j+1)%5],V[j]);h=diff(V[(j-1)%5],V[j]);el=iv.sqrt(dot(e,e))
        ux,uy=e[0]/el,e[1]/el
        ang=iv.atan2(cross(e,h),dot(e,h))
        if low(ang)<=0 or high(ang)>=math.pi: raise ValueError('not convex frame')
        nus=[(m+1)*P/ang for m in range(data['K'])]
        gs=[1/iv.gamma(nu+1) for nu in nus]
        frames.append((V[j],ux,uy,nus,gs))
    return V,k,coeffs,frames

def evaluate(point,edge,setup_data,order=2):
    V,k,coeffs,frames=setup_data
    total=iv.mpf(0);dtot=iv.mpf(0);ddtot=iv.mpf(0); dddtot=iv.mpf(0)
    for j in range(5):
        if edge is not None and (j==edge or j==(edge+1)%5):continue
        o,c,s,nus,gs=frames[j]
        w=diff(point,o);wx=c*w[0]+s*w[1];wy=-s*w[0]+c*w[1]
        rr=iv.sqrt(wx*wx+wy*wy);theta=iv.atan2(wy,wx)
        z=k*rr
        if edge is not None:
            E=diff(V[(edge+1)%5],V[edge]);e=E
            rp=dot(w,e)/rr; thp=cross(w,e)/(rr*rr)
            rpp=(dot(e,e)-rp*rp)/rr;thpp=-2*rp*thp/rr
            zp=k*rp;zpp=k*rpp
            if order>=3:
                rppp=-3*rp*rpp/rr
                thppp=thp*(-2*rpp/rr+6*rp*rp/(rr*rr))
                zppp=k*rppp
        for m in range(len(nus)):
            nu=nus[m];c0=coeffs[j*len(nus)+m]
            J,Jd,Jdd=bessel_series(nu,z,gs[m])
            S=iv.sin(nu*theta);C=iv.cos(nu*theta)
            total+=c0*J*S
            if edge is not None:
                dtot+=c0*(Jd*zp*S+J*C*nu*thp)
                if order>=2:
                    ddtot+=c0*((Jdd*zp*zp+Jd*zpp)*S + 2*Jd*zp*C*nu*thp+J*(-S*(nu*thp)**2+C*nu*thpp))
                if order>=3:
                    Jddd= -Jdd/z+Jd/(z*z) - (1 - nu*nu/(z*z))*Jd - 2*nu*nu*J/(z*z*z)
                    A1=Jd*zp; A2=Jdd*zp*zp+Jd*zpp
                    A3=Jddd*zp*zp*zp+3*Jdd*zp*zpp+Jd*zppp
                    b1=nu*thp*C
                    b2=nu*thpp*C - (nu*thp)**2*S
                    b3=(nu*thppp-(nu*thp)**3)*C-3*(nu*thp)*(nu*thpp)*S
                    dddtot+=c0*(A3*S+3*A2*b1+3*A1*b2+J*b3)
    return total,dtot,ddtot,dddtot

def run(inp,segments=100,edges=None,startseg=0,stopseg=None):
    start=time.time();data=json.load(open(inp));sd=setup(data);V,k,c,fr=sd
    # incenter disk radius: centroid to every edge supporting line, positive half planes separately
    x0=[sum(x[j] for x in V)/5 for j in (0,1)]
    radii=[]
    for j in range(5):
        a=V[j];b=V[(j+1)%5];e=diff(b,a)
        r=cross(e,diff(x0,a))/iv.sqrt(dot(e,e))
        if low(r)<=0:raise ValueError('centroid not inside')
        radii.append(low(r))
    R=iv.mpf(str(min(radii)*0.99))
    u0=evaluate(x0,None,sd)[0]
    # integer-order J1 via same proven series; nu=1, 1/gamma(2)=1
    J1=bessel_series(iv.mpf(1),k*R,iv.mpf(1),N=38)[0]
    intB=u0*2*P*R*J1/k
    # division by pi*R^2, square may be interval crossing 0; positive lower defined rigorously
    l2=(intB*intB/(P*R*R))
    l2_lb=max(0,low(l2))
    print('disk R',R,'u0',u0,'L2 >=',l2_lb,flush=True)
    supremum=0
    for edge in (edges if edges is not None else range(5)):
        a=V[edge];e=diff(V[(edge+1)%5],a)
        edge_sup=0
        for si in range(startseg,stopseg if stopseg is not None else segments):
            if si % 100 == 0: print('progress',edge,si,flush=True)
            tm=iv.mpf(2*si+1)/iv.mpf(2*segments)
            mid=add(a,scal(tm,e))
            # midpoint value and t derivative; second derivative evaluated on actual interval of t, rigorous
            f,df,ddmid,_=evaluate(mid,edge,sd,order=2)
            qlo=iv.mpf(si)/iv.mpf(segments);qhi=iv.mpf(si+1)/iv.mpf(segments)
            tt=iv.mpf([qlo.a,qhi.b])
            xx=add(a,scal(tt,e))
            _,_,dd,ddd=evaluate(xx,edge,sd,order=3)
            h=math.nextafter(1/(2*segments),math.inf)
            ub=abmax(f)+abmax(df)*h+abmax(ddmid)*h*h/2+abmax(ddd)*h*h*h/6
            edge_sup=max(edge_sup,ub)
        print('edge',edge,'sup',edge_sup,'time',round(time.time()-start,1),flush=True)
        supremum=max(supremum,edge_sup)
    area=iv.mpf(0)
    for j in range(5):area+=cross(V[j],V[(j+1)%5])/2
    mu=iv.sqrt(area)*iv.mpf([str(supremum),str(math.nextafter(supremum,math.inf))])/iv.sqrt(iv.mpf(str(l2_lb)))
    lam=Q(data['lam_tilde']);F=(lam/(1+mu))*area
    F_lower=low(F);hks=high(lam/(1-mu)*area) < 36.33682907107446
    print('RESULT',{'file':inp,'N':segments,'sup':supremum,'L2_lower':l2_lb,'mu_upper':high(mu),'F_lower':F_lower,'producer_F_lower':data.get('producer_F_lo_v6',data.get('producer_F_lo',0)),'gap':F_lower-float(data.get('producer_F_lo_v6',data.get('producer_F_lo',0))),'HKS':hks,'elapsed':time.time()-start},flush=True)
    return F_lower
if __name__=='__main__':
    run(sys.argv[1],int(sys.argv[2]) if len(sys.argv)>2 else 100, [int(sys.argv[3])] if len(sys.argv)>3 else None, int(sys.argv[4]) if len(sys.argv)>4 else 0, int(sys.argv[5]) if len(sys.argv)>5 else None)
