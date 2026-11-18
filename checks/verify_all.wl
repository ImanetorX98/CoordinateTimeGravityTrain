(* Symbolic verification of the analytic results of the paper, in Wolfram Language.
   Same content as checks/check_derivations.py (SymPy), independently implemented.
   Run:  wolframscript -file checks/verify_all.wl     (or open it in Mathematica)
   Output: a label and True/False for each check, then a summary line. *)
ClearAll["Global`*"];
$Assumptions = 0 < q < x < 1 && 0 < mu < 8/9;
passed = 0; failed = 0;
test[label_, value_] := Module[{ok = TrueQ[value]},
 Print[label]; Print[ok]; If[ok, passed++, failed++]];
eq[expr_] := TrueQ[FullSimplify[expr == 0, $Assumptions]];
tr[expr_, var_, order_] := Normal[Series[expr, {var, 0, order}]];

(* Sec. VI A-B: W must be retained through mu^3, despite n being required through mu^2. *)
ax = (3 Sqrt[1-mu]-Sqrt[1-mu x^2])^2/4;
ar = 1-mu;
aclaim = 1-mu (3-x^2)/2+3 mu^2 (1-x^2)^2/16;
wseries = Expand[tr[ax (ar-ax),mu,3]];
vn = Cancel[2 wseries/(mu (1-x^2))];
n1 = 7 (1-x^2)/16; n2 = (1-x^2) (227-3 x^2)/512;
normn = tr[Sqrt[(1-mu)/vn],mu,2];
test["Sec. VI A-B A through mu^2",eq[tr[ax-aclaim,mu,2]]];
test["Sec. VI A-B surface gap through mu^2",eq[tr[ar-ax-mu (1-x^2)/2 (1-3 mu (1-x^2)/8),mu,2]]];
test["Sec. VI A-B normalized n through mu^2",eq[normn-1-mu n1-mu^2 n2]];

(* Sec. VI C: exact series of squared normalized integrand; all branches positive. *)
mx = Expand[tr[(1+mu n1+mu^2 n2)^2/(1-x^2),mu,2]];
mq = mx /. x->q;
den = Expand[x^2 mx-q^2 mq];
d0 = Factor[Coefficient[den,mu,0]];
d1 = Factor[Coefficient[den,mu,1]/d0]; d2 = Factor[Coefficient[den,mu,2]/d0];
ratsq = Expand[tr[(1+mu x^2+mu^2 x^4) mq (x^2-q^2)/((1-x^2) d0) (1-mu d1+mu^2 (d1^2-d2)),mu,2]];
r1 = Factor[Coefficient[ratsq,mu,1]/2];
r2 = Factor[(Coefficient[ratsq,mu,2]-r1^2)/2];
i0=q Sqrt[1-x^2]/(x Sqrt[x^2-q^2]);
i1=q (15-7 q^2) x Sqrt[1-x^2]/(16 Sqrt[x^2-q^2]);
i2=3 q x Sqrt[1-x^2] ((49 q^4-118 q^2+133) x^2-48 q^4+4 q^2+44)/(512 Sqrt[x^2-q^2]);
test["Sec. VI C I0 from exact squared integrand",eq[Coefficient[ratsq,mu,0]-1]];
test["Sec. VI C I1/I0",eq[r1-i1/i0]];
test["Sec. VI C I2/I0",eq[r2-i2/i0]];

(* Sec. VI D-E: actual derivatives, zero lower endpoints, upper endpoints and reversion. *)
uangle=ArcSin[Sqrt[(x^2-q^2)/(1-q^2)]];
sroot=Sqrt[(x^2-q^2) (1-x^2)];
pp=-144 q^4+12 q^2+132; pq=147 q^4-354 q^2+399;
phi1=q (15-7 q^2)/32 (sroot+(1-q^2) uangle);
phi2=q/512 (pp/2 (sroot+(1-q^2) uangle)+pq/8 ((2 x^2+3 q^2-1) sroot+(1-q^2) (1+3 q^2) uangle));
test["Sec. VI D-E derivative U",eq[D[uangle,x]-x/sroot]];
test["Sec. VI D-E dphi1/dx = I1",eq[D[phi1,x]-i1]];
test["Sec. VI D-E dphi2/dx = I2",eq[D[phi2,x]-i2]];
test["Sec. VI D-E phi1(q)=phi2(q)=0",eq[phi1/.x->q] && eq[phi2/.x->q]];
f1=Pi q (1-q^2) (15-7 q^2)/64;
f2=3 Pi q (1-q^2) (147 q^6-497 q^4+297 q^2+309)/8192;
test["Sec. VI D-E integrated NLO endpoint",eq[(phi1/.x->1)-f1]];
test["Sec. VI D-E integrated NNLO endpoint",eq[(phi2/.x->1)-f2]];
test["Sec. VI D-E Newtonian primitive after u=U",eq[q/(q^2 Cos[u]^2+Sin[u]^2)-q-q (1-q^2) Cos[u]^2/(q^2+(1-q^2) Sin[u]^2)]];
q1=q (1-q^2) (15-7 q^2)/32;
pol=1827-3489 z+2457 z^2-539 z^3;
q2=q (1-q^2) (pol/.z->q^2)/4096;
test["Sec. VI D-E inversion q1",eq[-Pi q1/2+f1]];
test["Sec. VI D-E inversion q2",eq[-Pi q2/2+D[f1,q] q1+f2]];
bern=1827 (1-z)^3+1992 z (1-z)^2+960 z^2 (1-z)+256 z^3;
test["Sec. VI D-E positive polynomial certificate",eq[pol-bern] && And@@Thread[{15,8,1827,1992,960,256}>0]];

(* Sec. VII A: regular and weighted central limits. *)
ab=ax/.mu->8/9-eps;
test["Sec. VII A sqrt(1-mu) at Buchdahl",eq[tr[Sqrt[1/9+eps],eps,1]-1/3-3 eps/2]];
test["Sec. VII A central A",eq[tr[ab/.x->0,eps,2]-81 eps^2/16]];
test["Sec. VII A central n leading coefficient",Limit[eps^2 (1/9+eps)/((ab/.x->0) (1/9+eps-(ab/.x->0))),eps->0]==16/81];
test["Sec. VII A exact-in-mu small-x coefficient",eq[tr[ax,x,2]-(ax/.x->0)-mu (3 Sqrt[1-mu]-1) x^2/4]];
test["Sec. VII A uniform central approximation",eq[tr[ab/.{eps->tt ee,x->Sqrt[tt xx]},tt,2]-tt^2 (9 ee+8 xx/9)^2/16]];

(* Sec. VII B, circular orbits: eliminate yy at a double zero of the derivative polynomial. *)
dpoly=(3 aa-yy)^2 (yy-aa) (5 aa-yy);
fpoly=-15 aa^3 yy+9 aa^2 yy^2+14 aa^2+3 aa yy^3-12 aa yy-yy^4+2 yy^2;
critical=3645 z^3+3888 z^2-1080 z+64;
test["Sec. VII B (circular orbits) G derivative reduction",eq[(1-yy^2) D[dpoly,yy]+2 yy dpoly-2 (3 aa-yy) fpoly]];
test["Sec. VII B (circular orbits) resultant",eq[Resultant[fpoly,D[fpoly,yy],yy]+16 aa^2 (aa^2-1) (25 aa^2-1) (critical/.z->aa^2)]];
ycrit=aa (2025 aa^6+3987 aa^4-474 aa^2+48)/(3240 aa^6+2664 aa^4-334 aa^2+16);
test["Sec. VII B (circular orbits) common root modulo cubic",And@@(PolynomialRemainder[Numerator[Together[#/.yy->ycrit]],critical/.z->aa^2,aa]===0 & /@ {fpoly,D[fpoly,yy]})];
test["Sec. VII B (circular orbits) unique physical critical z",CountRoots[critical,{z,1/9,1}]==1];
(* Apply the G'=0 rule to the whole completed identity; retain G''. *)
stability=D[ll^2/(ar0 gg[x]),{x,2}]+ll^2 gg''[x]/(ar0 gg[x]^2);
test["Sec. VII B (circular orbits) stability: max stable, min unstable",eq[stability/.Derivative[1][gg][x]->0]];
test["Sec. VII B (circular orbits) logarithmic bottleneck primitive",FullSimplify[D[ArcSinh[x/Sqrt[dd/bb]]/Sqrt[bb],x]==1/Sqrt[dd+bb x^2],dd>0 && bb>0 && x>0]];
test["Sec. VII B (circular orbits) critical bottleneck exponent",FullSimplify[dd/(dd^3)^(1/2)==dd^(-1/2),dd>0]];

(* Sec. VII B, diameter: exact regularization, Jacobi equation, cubic, local optimality. *)
vreg=mm (3 aa-yy)^2 (5 aa-yy)/(16 ww (yy+aa));
hreg=Sqrt[mm] (3 aa-yy) Sqrt[(5 aa-yy)/(yy+aa)]/(4 yy Sqrt[ww]);
du[expr_] := -mm/(2 yy) D[expr,yy];
loghu=mm/(2 yy) (1/(3 aa-yy)+1/(2 (5 aa-yy))+1/(2 (yy+aa))+1/yy);
test["Sec. VII B (diameter) surface factorization",eq[dpoly/16-ww vreg (yy^2-aa^2)/mm]];
test["Sec. VII B (diameter) regular integrand logarithmic derivative",FullSimplify[du[hreg]/hreg-loghu==0,mm>0 && ww>0 && 1/3<aa<yy<1 && yy<3 aa]];
jac=D[ff[x] vv[x],{x,2}]-ff''[x] vv[x];
test["Sec. VII B (diameter) Jacobi reduction of order",eq[jac/.{vv''[x]->-2 ff'[x]/ff[x]^3,vv'[x]->1/ff[x]^2}]];
test["Sec. VII B (diameter) cubic finite-part integration by parts",eq[D[ff[x^2],{x,2}]/x^2-D[ff[x^2],x]/x^3-4 ff''[x^2]]];
test["Sec. VII B (diameter) regular cubic integrand",FullSimplify[4/3 Sqrt[1-uu] D[(1-uu)^(3/2) mf[uu],{uu,2}]-(mf[uu]-4 (1-uu) mf'[uu]+4/3 (1-uu)^2 mf''[uu])==0,0<uu<1]];
tdiff=nn0 (kk qr^2+3 lc qr^4/4);
test["Sec. VII B (diameter) faster bifurcating branch when cubic negative",eq[(tdiff/.{qr^4->4 kk^2/lc^2,qr^2->-2 kk/lc})-nn0 kk^2/lc]];

(* Numerical values are evaluations of exact conditions, not fitted coefficients. *)
zc=z/.FindRoot[critical==0,{z,14/100,15/100},WorkingPrecision->65,AccuracyGoal->55,PrecisionGoal->55];
ac=Sqrt[zc]; yc=N[ycrit/.aa->ac,60]; muc=1-zc; xc=Sqrt[(1-yc^2)/muc];
Print["VALUE Sec. VII B mu_c = ",N[muc,48]]; Print["VALUE Sec. VII B x_c = ",N[xc,48]];
test["Sec. VII B (circular orbits) NUMERIC residual and physical domain",Abs[fpoly/.{aa->ac,yy->yc}]<10^-48 && Abs[D[fpoly,yy]/.{aa->ac,yy->yc}]<10^-48 && 1/3<ac<yc<1];

(* h=cos(theta) H, u=sin²(theta). K=integral[-H+2(1-u)H_u]. *)
wcentral[mm_] := (3 Sqrt[1-mm]-1)^2 (1-Sqrt[1-mm]) (5 Sqrt[1-mm]-1)/16;
kexpr=hreg (-1+2 (1-uu) loghu)/.{aa->Sqrt[1-mm],ww->wcentral[mm],yy->Sqrt[1-mm uu]};
kern[theta_,mval_] = kexpr/.{uu->Sin[theta]^2,mm->mval};
kval[mval_?NumericQ,prec_:65] := NIntegrate[kern[theta,mval],{theta,0,Pi/4,Pi/2},WorkingPrecision->prec,AccuracyGoal->prec-15,PrecisionGoal->prec-15];
ms=mm/.FindRoot[kval[mm]==0,{mm,76/100,78/100},WorkingPrecision->65,AccuracyGoal->48,PrecisionGoal->48];
Print["VALUE Sec. VII B mu_s = ",N[ms,48]];
test["Sec. VII B (diameter) NUMERIC analytic condition",Abs[kval[ms]]<10^-45 && 7685/10000<ms<7686/10000];
(* Independent finite-part representation, with endpoint cancellation stabilized. *)
kfp[mval_?NumericQ] := Module[{aval=Sqrt[1-mval],w0=wcentral[mval],a0,a2,a4,w2,w4,h2,h4,fun},
 a0=(3 aval-1)^2/4; a2=mval (3 aval-1)/4; a4=3 aval mval^2/16;
 w2=a2 (aval^2-2 a0); w4=a4 (aval^2-2 a0)-a2^2;
 h2=(mval+w2/w0)/2; h4=(mval^2+mval w2/w0+w4/w0)/2-(mval+w2/w0)^2/8;
 fun[theta_?NumericQ] := If[Abs[Sin[theta]]<10^-12,(h2+h4 Sin[theta]^2) Cos[theta],
 ((Cos[theta] (hreg/.{aa->aval,ww->w0,mm->mval,yy->Sqrt[1-mval Sin[theta]^2]}))-1) Cos[theta]/Sin[theta]^2];
 NIntegrate[fun[theta],{theta,0,10^-12,Pi/4,Pi/2},WorkingPrecision->65,AccuracyGoal->40,PrecisionGoal->40]-1];
test["Sec. VII B (diameter) NUMERIC independent finite-part condition",Abs[kfp[ms]]<10^-35];
mexpr=mm^(3/2) (3 aa-yy)^3 (5 aa-yy)^(3/2)/(64 ww^(3/2) yy (yy+aa)^(3/2));
jexpr=mexpr-4 (1-uu) du[mexpr]+4/3 (1-uu)^2 du[du[mexpr]];
jkernel[theta_,mval_] = (jexpr/.{aa->Sqrt[1-mm],ww->wcentral[mm],yy->Sqrt[1-mm uu]})/.{uu->Sin[theta]^2,mm->mval};
js=NIntegrate[jkernel[theta,ms],{theta,0,Pi/4,Pi/2},WorkingPrecision->65,AccuracyGoal->45];
kder[theta_,mval_] = D[kern[theta,mm],mm]/.mm->mval;
ks=NIntegrate[kder[theta,ms],{theta,0,Pi/4,Pi/2},WorkingPrecision->65,AccuracyGoal->45];
Print["VALUE Sec. VII B K_prime = ",N[ks,32]]; Print["VALUE Sec. VII B cubic = ",N[js,32]];
test["Sec. VII B (diameter) NUMERIC transverse crossing and supercritical branch",ks>0 && js<0];
ms75=mm/.FindRoot[kval[mm,75]==0,{mm,76/100,78/100},WorkingPrecision->75,AccuracyGoal->58,PrecisionGoal->58];
test["Sec. VII B (diameter) NUMERIC precision stability",Abs[ms-ms75]<10^-45];

(* Sec. VI G, Earth example: rounded physical inputs quoted in the paper. *)
qearth=N[1-(522/100000)/Pi,50]; mearth=139/10^11; rearth=6371000; cspeed=299792458;
eq1=qearth (1-qearth^2) (15-7 qearth^2)/32;
shift=rearth mearth eq1; omega=Sqrt[mearth cspeed^2/(2 rearth^2)];
tearth=Pi/omega Sqrt[1-qearth^2];
Print["VALUE Sec. VI G q1 = ",N[eq1,18]," depth reduction (micrometres) = ",N[shift 10^6,18]," time(s) = ",N[tearth,18]];
test["Sec. VI G (Earth example) NUMERIC rounded Earth example",Abs[eq1-831/10^6]<2/10^6 && Abs[shift 10^6-74/10]<1/10 && Abs[tearth-146]<1];
test["Sec. VI G (Earth example) omega identity",eq[(mu cc^2/(2 rr^2)/.mu->2 gm/(rr cc^2))-gm/rr^3]];
Print["SUMMARY: ",passed," PASS; ",failed," FAIL; 0 SKIP"];
If[failed>0,Exit[1]];
