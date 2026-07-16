# t-branch clock: ENTIRELY ON-CURVE (frame-drag cover era spurio)

## Correzione
Il paper (eq:clock-t vecchia) dava dt/dr = ρ_t/√R6 + c_β√(2Mr/(r²+a²)) con un
termine di **doppio ricoprimento** (frame-dragging cover). **È SBAGLIATO.**
Il clock del t-branch è **interamente sulla curva** genus-2 R6, come il τ-branch.

## Identità chiave (esatta, sympy)
    (r-2M)Q2 + DE(J(r-2M)+2Ma)² = E² r³ Δ
con DE=(E²-1)r+2M, e  E²-f = DE/r.  Quindi il termine di Randers
    n_t α = E² r³/(f √R6)   ESATTAMENTE (nessun ricoprimento).

## Forma chiusa corretta
    dt/dr = ρ_t/√R6,   ρ_t = [E² r³ - (2Ma/r) K]/f = P_3(r) + R_Δ(r)/Δ
- polo all'ergosfera r=2M **rimovibile** (numeratore si annulla lì)
- poli terza specie SOLO agli orizzonti (Δ=0)

**Parte seconda specie (polinomio):**
    P_3 = E²r³ + 2E²M r² + 4E²M² r + 8E²M³ - 2(E²-1)JMa
**Vettore clock b** = (8E²M³-2(E²-1)JMa, 4E²M², 2E²M, E², 0)
  (analogo t del τ-branch b=(0,0,-2M,1,0))
**Parte terza specie orizzonte** (numeratore grado 1):
    R_Δ = -2M(2E²JMar - E²Ja³ - 8E²M³r + 4E²M²a² + 2E²Ma²r + Ja³ - 2Ma²r)

## Verifica
- Identità G=E²r³Δ: esatta (sympy)
- dt/dr on-curve vs **flusso di Hamilton** (H del ramo t): match 1e-10 a r=4,6,9,11
  (segno = solo orientamento ingoing)
- Perché il vecchio "verificato 1e-14" col cover era sbagliato: c_β√(2Mr/(r²+a²))
  a r=6 vale ~0.85 (NON piccolo); se fosse nel dt/dr, sballerebbe di 0.85, ma la
  forma on-curve combacia col flusso a 1e-10. Radicali diversi (√R6 vs cover) non
  possono coincidere con c_β≠0.

## Conseguenza (semplificazione)
Il t-branch è ora **pienamente parallelo** al τ-branch: entrambi on-curve, δφ = ψ
assemblato coi W_kj di R6 + (t-branch) log-σ terza specie agli orizzonti. Niente
doppio ricoprimento. Il τ ha dipolo peso-2 all'∞; il t ha lettera terza specie
agli orizzonti.

Script: `reproduce_reductions.py` (identità + P_3 + b + R_Δ),
`kerr_tbranch_separatrix_weierstrass.py` (Q2, flusso di Hamilton).
Paper: eq:clock-t corretta (main + PRD), tab:branch-adiab, assemblaggio ψ.

## Chiusura del cerchio: assemblaggio psi t-branch (end-to-end)
Clock t = eta2 + eta3:  eta2 = sum_j b_j U_j (2a specie), eta3 = int R_Delta/(Delta sqrt R6) (3a specie orizzonte).
- Assemblaggio eq:psi-split (parte 2a specie): I_direct = I_assembled = 10.31674453,
  diff **8.85e-13** (riorganizzazione W_kj di R6 col vettore b^t).
- Clock completo: int dEF*t = I_assembled + int dEF*eta3 = 10.41916368, diff **8.85e-13**.
- pezzo orizzonte int dEF*eta3 = 0.1024 (log-sigma peso-2 a r_pm).
Orbita a J=2.5 = scattering (nessun turning), integrata fino a rp+0.6 (presso orizzonte).
Script: `kerr_tbranch_psi_assembly.py`. Paper: assemblaggio ψ + tab:script-map.
