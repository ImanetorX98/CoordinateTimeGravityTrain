# Progresso: forme chiuse esplicite delle brachistocrone adiabatiche

Stato al: sessione in corso. Obiettivo: **forma chiusa esplicita in FUNZIONI SPECIALI**
(anche non elementari) delle correzioni adiabatiche al primo ordine δφ, testabile
contro la numerica. Metodo: **tutto analitico, coefficienti da residui/algebra, mai
fit; ogni passo verificato numericamente prima di costruirci sopra.**

---

## 1. Mappa 2×2 e stato

| | Vaidya (a=0) | Thakurta-Kerr (a≠0) |
|---|---|---|
| **J generico** (genus 2) | δφ assemblato (v,τ) ✓ rappr. integrale, testato 1e-13/1e-15. Forma in funz. speciali iperellittiche = DA FARE | δφ assemblato (t,τ) ✓ testato 5e-14/9e-13. Idem |
| **Separatrice** (genus 1) | **← QUI.** Settore abeliano/2ª specie ESPLICITO in σ,ζ,℘ ✓. Weight-2 = polilog ellittico (frontiera) | non iniziato (stessa macchina) |

**Brachistocrone adiabatiche al 1° ordine — tutte fatte e verificate (rappr. integrale):**
- TK τ: 5e-14 | TK t (on-curve, cover spurio corretto): 9e-13
- Vaidya v (tempo avanzato): 2e-13 | Vaidya τ (tempo proprio): 9.77e-15

---

## 2. Separatrice Vaidya: dati della curva (M=1, E=1.4, J variabile; separatrice a Jc)

- **Jc = 7.02662374** (radice doppia bracket), **r_d = −3.3637111** (radice doppia)
- Curva ellittica **E: w²=Q4(r)**, radici **{−2.0833, 0, 2, 8.7274}**, a4=0.96
- **k²=0.60672, τ=0.9059733802550 i** (Legendre = Sage, 15 cifre)
- Semiperiodi: om1=0.66913 (reale), om_im=0.60621
- z_∞ (immagine ∞) = 0.23663; z_d (immagine r_d) = 0.46104; z_h (2-torsione orizzonte) = 0.60621 i
- c_r = 0.13406; ρ = 1/√Q4(r_d) = 0.06107; C0 = 0.33327
- Turning fisico sulla separatrice = e4 = 8.7274 → orbita r∈(8.73, 12], z reale
- σ,ζ,℘ da θ₁ (mpmath.jtheta), ℘ ESATTO da θ₁''; reticolo auto-consistente

---

## 3. FATTO: settore abeliano / 1ª-2ª-3ª specie (ESPLICITO, verificato)

Tutte funzioni speciali esplicite (Weierstrass σ,ζ,℘ + algebrico), NO integrali in
quadratura. Coefficienti analitici (residui + ricorsione Hermite). Verifica 1e-9…1e-14.

**Formule (LS ≡ ln[σ(z−z_∞)/σ(z+z_∞)]):**
- U₀ = ρ[lnσ(z−z_d) − lnσ(z+z_d)] + C0·z      (3ª specie a r_d) — 1e-14
- V₁ = c_r·z − (1/√a4)·LS                       (3ª specie ∞)
- V₂ = c_r²z − (2c_r/√a4)LS + (1/a4)[−2ζ(2z_∞)LS − ζ(z−z_∞)−ζ(z+z_∞) + C·z]  (2ª specie)
- V₃,V₄,V₅ : **ricorsione di Hermite** da d/dr(r^k√Q4):
    (2k+4)a4 V_{k+3}+(2k+3)b3 V_{k+2}+(2k+2)b2 V_{k+1}+(2k+1)b1 V_k+2k b0 V_{k-1}=2 r^k√Q4
  con [a4,b3,b2,b1,b0] = coeff di Q4 (np.poly(radici)*a4)
- **U_k = Σ_{i=0}^{k-1} r_d^{k-1-i} V_i + r_d^k U₀**   (dallo split r^k/(r−r_d))
- **Π_h = β ζ(z−z_h) + γ z**,  β=−4/Q4'(2m)=0.07584, γ=−0.198751  (2ª specie orizzonte) — 1e-14

**Script:** `VaidyaMetric/vaidya_separatrix_explicit_Uk.py` (U₀..U₅, verif 1e-9),
`VaidyaMetric/vaidya_separatrix_Pih.py` (Π_h, 1e-14),
`VaidyaMetric/vaidya_ell_dilog_match.py` (U₀ σ,ζ + toolkit).

---

## 4. FRONTIERA: settore weight-2 (polilog ellittico) — NON ancora chiuso

Struttura (fatta l'algebra): W_kj = 2∫U_k dU_j − U_k U_j; ogni termine si riduce a
σ,ζ,℘ TRANNE **∫lnσ(z−a) ζ(z−b) dz = DILOGARITMO ELLITTICO** (Brown-Levin Γ̃ / Zagier D^E).
Analogamente D_k = U_k ln(r−2m) − G_k, G_k=∫ln(r−2m) r^k/√S dr = polilog ellittico.

**Cosa ho implementato (framework, verificato parzialmente):**
- kernel di Kronecker g^(1)(ξ)=θ₁'/θ₁·π + termine Im (single-valued), dispari ✓
- Zagier D^E(z)=Σ_{n≥0}D(qⁿζ)−Σ_{n≥1}D(qⁿ/ζ), e Li₂^ell olomorfo
- coordinate normalizzate ξ_h=τ/2=0.453i, ξ_d=0.3445

**Cosa NON chiude:** il **regulator preciso** — G₀ (reale, ~−0.001…−0.006) non ha
relazione lineare semplice con D^E/Li₂^ell (che sono ordine 2-7). La combinazione
esatta (poli a z_h, z_d E ∞, + correzioni single-valued) richiede la **macchina
completa di Brown-Levin**. È la **frontiera research 2024** (Broedel, Zerbini,
Schottky-Kronecker). NON truccato: il check numerico dice che non torna.

**Script tentativi (framework):** codice inline nelle ultime run (non salvato come file
definitivo — il regulator non chiude ancora).

---

## 5. Findings / correzioni importanti (verificati)

- **L_2m è SECONDA specie**, NON Fay: 𝒜^m(2m)=−19.6≠0 → integrando ~(r−2m)^{−3/2},
  polo doppio senza residuo → Kleinian ζ. Il genuino 3ª specie orizzonte è la lettera
  dr/(r−2m)=tortoise ln (elementare). **Fay serve in TK (lettera ∞, iperellittica
  r²/√S, non elementare), NON in Vaidya (lettera orizzonte razionale → ln elementare).**
  Paper corretto.
- **t-branch TK: clock ON-CURVE**, il "cover" di frame-dragging era SPURIO. Identità
  (r−2M)Q2+DE(J(r−2M)+2Ma)²=E²r³Δ ⇒ n_t α=E²r³/(f√R6). Corretto + verificato 1e-10.
- **δφ|_sep NON si ottiene mettendo Jc nella riduzione genus-2** (SINGOLARE lì, de Rham
  degenera). Va fatto come limite J→Jc o riduzione ellittica diretta. (L_2m diretto a Jc
  dava ~1e10, assurdo.)
- 𝓘_poly reso esplicito: polinomio + Σ_ρ res·log(r−ρ) sui 6 poli di S, verif 1e-15.

---

## 6. Software per chiudere il weight-2 (in corso di installazione)

Il regulator/polilog ellittico esplicito richiede strumenti dedicati:
- **GiNaC** (C++, `brew install ginac` → 1.8.10 + dipendenza CLN) — INSTALLARE
- **eMPL** (elliptic MPL con argomenti arbitrari): arXiv **2602.09956**, in GiNaC.
  Codice negli ancillary files del paper o repo linkato — DA RECUPERARE
- (genus-2 futuro) **Schottky-Kronecker**: arXiv 2406.10051 (Broedel-Zerbini), ancillary
- PolyLogTools (gitlab.com/pltteam/plt) è genus-0, NON serve al dilog ellittico
Ambiente: clang++ ✓, brew ✓, pkg-config ✓. GiNaC/CLN da installare.

---

## 7. DA DOVE RIPRENDERE

1. **Installare GiNaC** (`brew install ginac`) + recuperare il codice **eMPL** (arXiv:2602.09956).
2. Usare eMPL per valutare/chiudere il **dilog ellittico** ∫lnσ ζ (il regulator) — l'unico
   pezzo mancante della separatrice. Verificare vs G₀ diretto.
3. Assemblare **δφ|_sep = [σ,ζ,℘ espliciti (§3)] + [polilog ellittico (§4)]** come formula
   unica, via **limite J→Jc** (non riduzione a Jc, §5). Testare vs numerica.
4. Trasportare a **J generico (genus 2)** con Schottky-Kronecker, e a **TK** (stessa macchina).
5. Scrivere nel paper la sezione separatrice tutta-esplicita con citazioni (Fay, Baker,
   Buchstaber-Enolskii-Leykin, Brown-Levin, Beilinson-Levin, Zagier, Schottky-Kronecker).

**In sintesi:** settore abeliano ESPLICITO e verificato (§3); manca solo il **dilog
ellittico** del weight-2 (§4), che è frontiera e richiede eMPL/GiNaC. δφ|_sep va
assemblato via limite J→Jc.
