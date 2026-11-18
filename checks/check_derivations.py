#!/usr/bin/env python3
"""Symbolic verification of the analytic results of the paper.

Purpose  Check every analytic formula of Secs. VI-VII by exact algebra (SymPy),
         plus high-precision evaluation (mpmath) of the two thresholds.
Paper    Sec. VI A-B   expansions of A, B and of the optical index n
         Sec. VI C     coefficients I0, I1, I2 of the Snell integrand (Eqs. 47-49)
         Sec. VI D-E   closed-form integrals, Eq. (54), and the inversion q1, q2 (Eqs. 57-59),
                       including a proof that q1 and q2 are positive
         Sec. VI G     the Earth example
         Sec. VII A    asymptotics near the Buchdahl bound
         Sec. VII B    onset of circular optical orbits, mu_c (Eq. 70), and loss of
                       local minimality of the diameter, mu_s (Eqs. 71-73)
Inputs   none.
Outputs  "PASS | ..." / "FAIL | ..." lines, "VALUE | ..." lines with the thresholds,
         and a summary; exit status 0 if all pass. Checks marked NUMERIC evaluate
         analytically derived conditions to many digits; nothing is fitted.
Runtime  about 10 s.
"""

import sys
import sympy as s
import mpmath as mp

results = []
def check(label, fn):
    try:
        ok = bool(fn())
        detail = ''
    except Exception as exc:
        ok, detail = False, f' ({type(exc).__name__}: {exc})'
    results.append((label, ok))
    print(f'{"PASS" if ok else "FAIL"} | {label}{detail}', flush=True)

def zero(e):
    return s.factor(s.cancel(s.trigsimp(e))) == 0

def trunc(e, v, n):
    return s.series(e, v, 0, n).removeO().expand()

x, q, mu, e, z = s.symbols('x q mu epsilon z', positive=True)
a, y = s.symbols('a y', positive=True)
A = (3*s.sqrt(1-mu)-s.sqrt(1-mu*x*x))**2/4
AR = 1-mu
Ac = 1-mu*(3-x*x)/2+3*mu**2*(1-x*x)**2/16
# W is needed through mu^3 to obtain the normalized index through mu^2.
W = trunc(A*(AR-A), mu, 4)
vn = s.cancel(2*W/(mu*(1-x*x)))
n1 = 7*(1-x*x)/16
n2 = (1-x*x)*(227-3*x*x)/512
N = trunc(s.sqrt((1-mu)/vn), mu, 3)
check('Sec. VI A-B A through mu^2', lambda: zero(trunc(A-Ac,mu,3)))
check('Sec. VI A-B surface gap through mu^2', lambda: zero(trunc(AR-A-mu*(1-x*x)/2*(1-3*mu*(1-x*x)/8),mu,3)))
check('Sec. VI A-B normalized n through mu^2 (uses W through mu^3)', lambda: zero(N-1-n1*mu-n2*mu**2))

# n^2=(2/mu) M(x). Use series ring arithmetic for the exact squared ratio.
M = trunc((1+mu*n1+mu**2*n2)**2/(1-x*x),mu,3)
Mq = M.subs(x,q)
# I^2/I0^2 = (1/(1-mu*x^2))*Mq*(x^2-q^2)/((1-x^2)*(x^2*M-q^2*Mq)).
den = s.expand(x*x*M-q*q*Mq)
d0 = s.factor(den.coeff(mu,0))
d1 = s.factor(den.coeff(mu,1)/d0)
d2 = s.factor(den.coeff(mu,2)/d0)
ratio2 = trunc((1+mu*x*x+mu**2*x**4)*Mq*(x*x-q*q)/((1-x*x)*d0)*(1-mu*d1+mu**2*(d1*d1-d2)),mu,3)
r1 = s.factor(ratio2.coeff(mu,1)/2)
r2 = s.factor((ratio2.coeff(mu,2)-r1*r1)/2)
I0 = q*s.sqrt(1-x*x)/(x*s.sqrt(x*x-q*q))
I1 = q*(15-7*q*q)/16*x*s.sqrt(1-x*x)/s.sqrt(x*x-q*q)
I2 = 3*q/512*x*s.sqrt(1-x*x)/s.sqrt(x*x-q*q)*((49*q**4-118*q*q+133)*x*x-48*q**4+4*q*q+44)
check('Sec. VI C I0 from exact squared integrand', lambda: zero(ratio2.coeff(mu,0)-1))
check('Sec. VI C I1/I0', lambda: zero(r1-I1/I0))
check('Sec. VI C I2/I0', lambda: zero(r2-I2/I0))

# On q<x<1 the radical branches are positive. Verify the two elementary
# primitives by their positive-branch derivatives, avoiding unsafe power rules.
S2=(x*x-q*q)*(1-x*x)
S=s.sqrt(S2)
U=s.asin(s.sqrt((x*x-q*q)/(1-q*q)))
Sp=s.diff(S2,x)/(2*S)
Up=x/S
# U'=x/S follows by squaring the positive radicals in the chain rule.
check('Sec. VI D-E derivative U (squared positive-branch identity)', lambda: zero(s.diff((x*x-q*q)/(1-q*q),x)**2/(4*((x*x-q*q)/(1-q*q))*(1-(x*x-q*q)/(1-q*q)))-Up**2))
P=-144*q**4+12*q*q+132
Q=147*q**4-354*q*q+399
# Derivatives written over S: I1 = q(15-7q²)x(1-x²)/(16 S).
dp1=q*(15-7*q*q)/32*(Sp+(1-q*q)*Up)
dp2=q/512*(P/2*(Sp+(1-q*q)*Up)+Q/8*(4*x*S+(2*x*x+3*q*q-1)*Sp+(1-q*q)*(1+3*q*q)*Up))
check('Sec. VI D-E dphi1/dx = I1', lambda: zero(dp1-q*(15-7*q*q)/16*x*(1-x*x)/S))
check('Sec. VI D-E dphi2/dx = I2', lambda: zero(dp2-3*q/512*x*(1-x*x)/S*((49*q**4-118*q*q+133)*x*x-48*q**4+4*q*q+44)))
phi1=q*(15-7*q*q)/32*(S+(1-q*q)*U)
phi2=q/512*(P/2*(S+(1-q*q)*U)+Q/8*((2*x*x+3*q*q-1)*S+(1-q*q)*(1+3*q*q)*U))
check('Sec. VI D-E phi1(q)=phi2(q)=0', lambda: phi1.subs(x,q)==0 and phi2.subs(x,q)==0)
f1=s.pi/64*q*(1-q*q)*(15-7*q*q)
f2=3*s.pi/8192*q*(1-q*q)*(147*q**6-497*q**4+297*q*q+309)
check('Sec. VI D-E integrated NLO endpoint', lambda: zero(phi1.subs(x,1)-f1))
check('Sec. VI D-E integrated NNLO endpoint', lambda: zero(phi2.subs(x,1)-f2))
# Newtonian primitive: arctan(sqrt(x²-q²)/(q sqrt(1-x²)))-q U.
# After u=U, I0 dx=q(1-q²)cos²u/(q²+(1-q²)sin²u) du.
u=s.symbols('u', real=True)
check('Sec. VI D-E Newtonian primitive derivative after u=U', lambda: zero(q/(q*q*s.cos(u)**2+s.sin(u)**2)-q-q*(1-q*q)*s.cos(u)**2/(q*q+(1-q*q)*s.sin(u)**2)))
q1=q*(1-q*q)*(15-7*q*q)/32
poly=1827-3489*z+2457*z*z-539*z**3
q2=q*(1-q*q)*poly.subs(z,q*q)/4096
check('Sec. VI D-E inversion q1', lambda: zero(-s.pi*q1/2+f1))
check('Sec. VI D-E inversion q2', lambda: zero(-s.pi*q2/2+s.diff(f1,q)*q1+f2))
bern=1827*(1-z)**3+1992*z*(1-z)**2+960*z*z*(1-z)+256*z**3
check('Sec. VI D-E positivity certificate for q1 and q2 on (0,1)', lambda: zero(poly-bern) and all(c>0 for c in [15,8,1827,1992,960,256]))

AB=A.subs(mu,s.Rational(8,9)-e)
check('Sec. VII A sqrt(1-mu) at Buchdahl', lambda: zero(trunc(s.sqrt(s.Rational(1,9)+e),e,2)-s.Rational(1,3)-3*e/2))
check('Sec. VII A central A = 81 epsilon^2/16', lambda: zero(trunc(AB.subs(x,0),e,3)-81*e*e/16))
check('Sec. VII A central n leading coefficient 4/9', lambda: s.limit(e**2*AR.subs(mu,s.Rational(8,9)-e)/(AB.subs(x,0)*(AR.subs(mu,s.Rational(8,9)-e)-AB.subs(x,0))),e,0)==s.Rational(16,81))
check('Sec. VII A exact-in-mu small-x coefficient', lambda: zero(trunc(A,x,3)-A.subs(x,0)-mu*(3*s.sqrt(1-mu)-1)*x*x/4))
# Weighted scaling epsilon=t E, x²=t X: leading uniform central square.
t,E,X=s.symbols('t E X',positive=True)
check('Sec. VII A uniform central approximation (weighted expansion)', lambda: zero(trunc(AB.subs({e:t*E,x:s.sqrt(t*X)}),t,3)-t*t*(9*E+s.Rational(8,9)*X)**2/16))

D=(3*a-y)**2*(y-a)*(5*a-y)
F=-15*a**3*y+9*a*a*y*y+14*a*a+3*a*y**3-12*a*y-y**4+2*y*y
H=(1-y*y)*s.diff(D,y)+2*y*D
crit=3645*z**3+3888*z*z-1080*z+64
check('Sec. VII B (circular orbits) G derivative reduction to F', lambda: zero(H-2*(3*a-y)*F))
res=s.resultant(F,s.diff(F,y),y)
check('Sec. VII B (circular orbits) exact resultant', lambda: zero(res+16*a*a*(a*a-1)*(25*a*a-1)*crit.subs(z,a*a)))
yc=a*(2025*a**6+3987*a**4-474*a*a+48)/(3240*a**6+2664*a**4-334*a*a+16)
check('Sec. VII B (circular orbits) critical common root F=Fy=0 modulo cubic', lambda: all(s.rem(s.together(f.subs(y,yc)).as_numer_denom()[0],crit.subs(z,a*a),a)==0 for f in [F,s.diff(F,y)]))
check('Sec. VII B (circular orbits) unique physical critical z in (1/9,1)', lambda: s.Poly(crit,z).count_roots(s.Rational(1,9),1)==1)
# Effective potential Veff=L²/(AR G); Veff''=-L² G''/(AR G²) at extrema.
g,L,ar=s.symbols('g L ar',positive=True)
gfun=s.Function('g')(x)
check('Sec. VII B (circular orbits) stability sign: max G stable, min G unstable', lambda: zero(s.diff(L**2/(ar*gfun),x,2).xreplace({s.diff(gfun,x):s.S.Zero})+L**2*s.diff(gfun,x,2)/(ar*gfun**2)))
d,b,eta=s.symbols('d b eta',positive=True)
check('Sec. VII B (circular orbits) logarithmic bottleneck primitive', lambda: zero(s.diff(s.asinh(x/s.sqrt(d/b))/s.sqrt(b),x)-1/s.sqrt(d+b*x*x)))
check('Sec. VII B (circular orbits) critical cubic bottleneck exponent -1/2', lambda: zero((d**3)**(-s.Rational(1,2))*d-d**(-s.Rational(1,2))))

# Sec. VII B (diameter): exact algebra for regularized finite-part integrals and Jacobi field.
# u=x²; V=W/[W0(1-u)], h=sqrt(1-u) H, H=sqrt(V)/y.
# W0 and mu are independent symbols during the y differentiations.
w=s.symbols('w',positive=True)
V=mu*(3*a-y)**2*(5*a-y)/(16*w*(y+a))
HH=s.sqrt(mu)*(3*a-y)*s.sqrt((5*a-y)/(y+a))/(4*y*s.sqrt(w))
Dy=lambda expr: -mu/(2*y)*s.diff(expr,y)
logHu=mu/(2*y)*(1/(3*a-y)+1/(2*(5*a-y))+1/(2*(y+a))+1/y)
check('Sec. VII B (diameter) surface factorization W=(1-x²)W0 V', lambda: zero(D/16-w*V*(y*y-a*a)/mu))
check('Sec. VII B (diameter) regular integrand H_u/H', lambda: zero(Dy(HH)/HH-logHu))
ff=s.Function('f')(x)
JJ=s.Function('j')(x)
# j=f integral ds/f² solves j''-(f''/f)j=0; check reduction of order.
vv=s.Function('v')(x)
check('Sec. VII B (diameter) Jacobi reduction of order', lambda: zero((s.diff(ff*vv,x,2)-s.diff(ff,x,2)*vv).subs(s.diff(vv,x,2),-2*s.diff(ff,x)/ff**3).subs(s.diff(vv,x),1/ff**2)))
fu=s.Function('F')(u)
check('Sec. VII B (diameter) cubic finite-part integration by parts identity', lambda: zero((s.diff(fu.subs(u,x*x),x,2)/x**2-s.diff(fu.subs(u,x*x),x)/x**3)-4*s.diff(fu,u,2).subs(u,x*x)))
MM=s.Function('M')(u)
check('Sec. VII B (diameter) regular cubic integrand after x=sin(theta)', lambda: zero(s.Rational(4,3)*s.sqrt(1-u)*s.diff((1-u)**s.Rational(3,2)*MM,u,2)-(MM-4*(1-u)*s.diff(MM,u)+s.Rational(4,3)*(1-u)**2*s.diff(MM,u,2))))
k,ell,n0=s.symbols('k ell n0',nonzero=True)
qq=s.symbols('qq')
tdiff=n0*(k*qq**2+s.Rational(3,4)*ell*qq**4)
check('Sec. VII B (diameter) bifurcating branch is faster when cubic < 0', lambda: zero(tdiff.subs(qq**2,-2*k/ell)-n0*k*k/ell))

# Numerical evaluation of analytic roots/definite integrals, never a curve fit.
mp.mp.dps=55
zc=mp.findroot(lambda zz:3645*zz**3+3888*zz**2-1080*zz+64,('.14','.15'))
ac=mp.sqrt(zc)
yy=ac*(2025*zc**3+3987*zc**2-474*zc+48)/(3240*zc**3+2664*zc**2-334*zc+16)
muc=1-zc
xc=mp.sqrt((1-yy**2)/muc)
print('VALUE | Sec. VII B mu_c =',mp.nstr(muc,48))
print('VALUE | Sec. VII B x_c  =',mp.nstr(xc,48))
fnum=s.lambdify((a,y),F,'mpmath')
fynum=s.lambdify((a,y),s.diff(F,y),'mpmath')
check('Sec. VII B (circular orbits) NUMERIC critical root residual and admissible domain', lambda: abs(fnum(ac,yy))<mp.mpf('1e-48') and abs(fynum(ac,yy))<mp.mpf('1e-48') and ac<yy<1 and ac>mp.mpf(1)/3)

def K(m):
    aa=mp.sqrt(1-m)
    w0=(3*aa-1)**2*(1-aa)*(5*aa-1)/16
    def integrand(theta):
        uu=mp.sin(theta)**2
        yy=mp.sqrt(1-m*uu)
        h=mp.sqrt(m)*(3*aa-yy)*mp.sqrt((5*aa-yy)/(yy+aa))/(4*yy*mp.sqrt(w0))
        lu=m/(2*yy)*(1/(3*aa-yy)+1/(2*(5*aa-yy))+1/(2*(yy+aa))+1/yy)
        return h*(-1+2*(1-uu)*lu)
    return mp.quad(integrand,[0,mp.pi/4,mp.pi/2])

mus=mp.findroot(K,('.76','.78'))
print('VALUE | Sec. VII B mu_s =',mp.nstr(mus,48))
check('Sec. VII B (diameter) NUMERIC exact analytic condition K(mu_s)=0', lambda: abs(K(mus))<mp.mpf('1e-48') and mp.mpf('.7685')<mus<mp.mpf('.7686'))
# Independently evaluate the original finite-part formula at two precisions.
def Kfp(m):
    aa=mp.sqrt(1-m); w0=(3*aa-1)**2*(1-aa)*(5*aa-1)/16
    # At small x use the exact convergent Taylor expansion through x^4;
    # this only stabilizes quadrature below 1e-12 (error O(1e-60)).
    av=(3*aa-1)**2/4
    a2=m*(3*aa-1)/4
    a4=3*aa*m*m/16
    w2=a2*(aa*aa-2*av)
    w4=a4*(aa*aa-2*av)-a2*a2
    h2=(m+w2/w0)/2
    h4=(m*m+m*w2/w0+w4/w0)/2-(m+w2/w0)**2/8
    def integrand(theta):
        xx=mp.sin(theta)
        if abs(xx)<mp.mpf('1e-12'):
            return (h2+h4*xx*xx)*mp.cos(theta)
        yy=mp.sqrt(1-m*xx*xx)
        # Rationalize y-a at the surface.
        hh=mp.cos(theta)*mp.sqrt(m)*(3*aa-yy)*mp.sqrt((5*aa-yy)/(yy+aa))/(4*yy*mp.sqrt(w0))
        return (hh-1)*mp.cos(theta)/(xx*xx)
    return mp.quad(integrand,[0,mp.mpf('1e-12'),mp.pi/4,mp.pi/2])-1

check('Sec. VII B (diameter) NUMERIC independent finite-part condition', lambda: abs(Kfp(mus))<mp.mpf('1e-35'))
Mexpr=mu**s.Rational(3,2)*(3*a-y)**3*(5*a-y)**s.Rational(3,2)/(64*w**s.Rational(3,2)*y*(y+a)**s.Rational(3,2))
mf=[s.lambdify((y,a,mu,w),expr,'mpmath') for expr in [Mexpr,Dy(Mexpr),Dy(Dy(Mexpr))]]
def cubic_J(m):
    aa=mp.sqrt(1-m); ww=(3*aa-1)**2*(1-aa)*(5*aa-1)/16
    def integrand(theta):
        yy=mp.sqrt(1-m*mp.sin(theta)**2); c2=mp.cos(theta)**2
        return mf[0](yy,aa,m,ww)-4*c2*mf[1](yy,aa,m,ww)+mp.mpf(4)/3*c2*c2*mf[2](yy,aa,m,ww)
    return mp.quad(integrand,[0,mp.pi/4,mp.pi/2])
js=cubic_J(mus)
kp=mp.diff(K,mus)
print('VALUE | Sec. VII B K_prime(mu_s) =',mp.nstr(kp,32))
print('VALUE | Sec. VII B cubic at mu_s =',mp.nstr(js,32))
check('Sec. VII B (diameter) NUMERIC transverse crossing and supercritical faster branch', lambda: kp>0 and js<0)
with mp.workdps(75):
    mus75=mp.findroot(K,('.76','.78'))
check('Sec. VII B (diameter) NUMERIC root stable at 55 and 75 decimal digits', lambda: abs(mus-mus75)<mp.mpf('1e-48'))

# Sec. VI G (Earth example). Constants are the rounded inputs quoted in the paper.
delta=mp.mpf('0.00522'); earth_mu=mp.mpf('1.39e-9'); R=mp.mpf('6371000'); c=mp.mpf('299792458')
q0=1-delta/mp.pi
earth_q1=q0*(1-q0*q0)*(15-7*q0*q0)/32
shift=R*earth_mu*earth_q1
omega=mp.sqrt(earth_mu*c*c/(2*R*R))
time=mp.pi/omega*mp.sqrt(1-q0*q0)
print('VALUE | Sec. VI G q1 =',mp.nstr(earth_q1,18),' depth reduction (micrometres) =',mp.nstr(shift*10**6,18),' time (s) =',mp.nstr(time,18))
check('Sec. VI G (Earth example) NUMERIC q1, depth reduction and time with rounded inputs', lambda: abs(earth_q1-mp.mpf('8.31e-4'))<mp.mpf('2e-6') and abs(shift*10**6-mp.mpf('7.4'))<mp.mpf('.1') and abs(time-146)<1)
gm,rr,cc=s.symbols('GM R c',positive=True)
check('Sec. VI G (Earth example) travel-time identity omega²=mu c²/(2R²)', lambda: zero((mu*cc**2/(2*rr**2)).subs(mu,2*gm/(rr*cc**2))-gm/rr**3))

passed=sum(ok for _,ok in results)
print(f'SUMMARY | {passed}/{len(results)} PASS; {len(results)-passed} FAIL; 0 SKIP')
sys.exit(0 if passed==len(results) else 1)
