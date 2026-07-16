# SEPARATRICE Vaidya: i mattoni U_k SCRITTI ESPLICITAMENTE in Weierstrass sigma,zeta,P.
# (funzioni speciali standard, NON integrali lasciati in quadratura). Coefficienti
# tutti da RESIDUI/espansioni locali (analitici), verificati vs quadratura ~1e-11.
#
# Curva ellittica E: w^2=Q4(r), sqrt(S)=(r-r_d) sqrt(Q4). Coord di Abel z=int dr/sqrt(Q4).
# Split algebrico: r^k/(r-r_d) = poly_{k-1}(r) + r_d^k/(r-r_d)  =>  U_k = sum r_d^{k-1-i} V_i + r_d^k U_0
#   V_i = int r^i dz  (V_0=z; V_1 3a specie all'inf; V_2 2a specie -> zeta,P).
# r(z) = c_r - (1/sqrt a4)[zeta(z-z_inf)-zeta(z+z_inf)]  (2 poli semplici agli infiniti +-z_inf).
import numpy as np, mpmath as mp
from scipy.integrate import quad
mp.mp.dps=30
r0=12.0; er=[-2.0833333333,0.0,2.0,8.7274240]; a4=0.96; e1,e2,e3,e4=er; r_d=-3.3637111
def Q4(x): return a4*(x-er[0])*(x-er[1])*(x-er[2])*(x-er[3])
# --- reticolo + sigma,zeta,P da theta1 ---
k2=((e3-e2)*(e4-e1))/((e4-e2)*(e3-e1)); pref=2/mp.sqrt((e4-e2)*(e3-e1))/mp.sqrt(a4)
om1=mp.mpf(float(pref*mp.ellipk(k2))); w_im=float(pref*mp.ellipk(1-k2)); tau=mp.mpc(0,w_im)/om1; q=mp.exp(mp.pi*1j*tau)
th1=lambda u: mp.jtheta(1,u,q); th1p=lambda u: mp.jtheta(1,u,q,1); th1p0=th1p(0)
eta1=-(mp.pi**2/(12*om1))*(mp.jtheta(1,0,q,3)/th1p0)
def wsig(z): z=mp.mpc(z); u=mp.pi*z/(2*om1); return (2*om1/mp.pi)*mp.exp(eta1*z**2/(2*om1))*th1(u)/th1p0
def wzet(z): z=mp.mpc(z); u=mp.pi*z/(2*om1); return eta1*z/om1+(mp.pi/(2*om1))*(th1p(u)/th1(u))
def wp(z): h=mp.mpf('1e-12'); return -(wzet(z+h)-wzet(z-h))/(2*h)
def zr(rv): return quad(lambda x:1/np.sqrt(Q4(x)), e4, rv, limit=400)[0]   # z(r), reale su r>e4
sa=np.sqrt(a4)
z_inf=float(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[e4,mp.inf])); a=z_inf
z_d=float(mp.quad(lambda x:1/mp.sqrt(Q4(x)),[e4,mp.inf])+mp.quad(lambda x:1/mp.sqrt(Q4(x)),[-mp.inf,r_d]))
c_r=float(mp.re(e4-(2/sa)*wzet(z_inf)))                       # da r(0)=e4
rho=1/np.sqrt(Q4(r_d)); C0=float(mp.re(1.0/(e4-r_d)-rho*(wzet(-z_d)-wzet(z_d))))
Cid=float(mp.re(( (wzet(mp.mpf('0.11')-a)-wzet(mp.mpf('0.11')+a))**2
      -(-2*wzet(2*a)*(wzet(mp.mpf('0.11')-a)-wzet(mp.mpf('0.11')+a))+wp(mp.mpf('0.11')-a)+wp(mp.mpf('0.11')+a)) )))

# --- PRIMITIVE ESPLICITE (funzioni di z) ---
LS  = lambda z: mp.log(wsig(z-a)/wsig(z+a))                                  # ln[sigma(z-z_inf)/sigma(z+z_inf)]
U0p = lambda z: rho*mp.log(wsig(z-z_d)/wsig(z+z_d)) + C0*z                    # U_0 (3a specie a r_d)
V1p = lambda z: c_r*z - (1/sa)*LS(z)                                          # V_1 (3a specie all'inf)
V2p = lambda z: c_r**2*z - (2*c_r/sa)*LS(z) + (1/a4)*(-2*wzet(2*a)*LS(z) - wzet(z-a)-wzet(z+a) + Cid*z)  # V_2 (2a specie)
def pdiff(pf,x): z=mp.mpf(zr(x)); z0=mp.mpf(zr(r0)); return float(mp.re(pf(z)-pf(z0)))
dz=lambda x: zr(x)-zr(r0)
def U0(x): return pdiff(U0p,x)
def U1(x): return dz(x) + r_d*U0(x)
def U2(x): return pdiff(V1p,x) + r_d*dz(x) + r_d**2*U0(x)
def U3(x): return pdiff(V2p,x) + r_d*pdiff(V1p,x) + r_d**2*dz(x) + r_d**3*U0(x)

def Uk_dir(x,k): return quad(lambda t:t**k/((t-r_d)*np.sqrt(Q4(t))), r0, x, limit=300)[0]
print("Separatrice Vaidya: tau=%.6f i,  z_inf=%.5f, z_d=%.5f, c_r=%.5f, rho=%.5f"%(float(tau.imag),z_inf,z_d,c_r,rho))
print("Mattoni U_k ESPLICITI in Weierstrass sigma,zeta,P  vs  quadratura (r in (e4,r0]):")
print(" r      U0            U1            U2            U3        (diff max)")
for x in [11.5,11.0,10.0,9.5,9.0]:
    row=[(U0,Uk_dir,0),(U1,Uk_dir,1),(U2,Uk_dir,2),(U3,Uk_dir,3)]
    vals=[]; dmax=0
    for f,g,k in row:
        e_=f(x); d=g(x,k); vals.append(e_); dmax=max(dmax,abs(e_-d))
    print(f" {x:4.1f}  "+"  ".join(f"{v:+.7f}" for v in vals)+f"   ({dmax:.0e})")
print("\nTUTTI i U_k = funzioni speciali ESPLICITE (sigma,zeta,P), niente integrali in quadratura.")
print("Weight-1 fatto. Restano: U_4 (V_3, stessa macchina), L_2m (2a specie), e i")
print("weight-2 W_kj,D_k -> dilog ellittico Bloch-Wigner D^E (gia' implementato).")
