# Admitted known forcings that remain distinguishable after a stiff–sloppy reduction

**Thesis #28.** Computational research, set out in Nile University B.Sc. chapter order for handoff.

**Depends on:** [Thesis #19](https://github.com/cloudynirvana/thesis-19-forcing-admission-gates-tip-ode) (forcing admission under evidence gates) and [Thesis #24](https://github.com/cloudynirvana/thesis-24-reduction-preserving-multichannel-id) (reduction-preserving multi-channel identifiability).

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**GitHub:** https://github.com/cloudynirvana  
**Date:** 21 September 2026

After a documented stiff–sloppy or MBAM-style reduction that preserves multi-channel ranks, which Thesis #19–admitted known forcings stay distinguishable from soft-prior leakage into Θ, and which admissions become artefacts of the unreduced coordinates?

The product forcing and the pool forcing stay distinguishable. The gain forcing is confounded with a coordinate write on both charts. The readout forcing is an artefact of the unreduced channel. The shared practical rank is 5 of 8 on the full vector and 5 of 6 on the reduced vector, so the reduction meets the rank condition before any forcing is scored. A one-coordinate leakage search on the factor interval [1/4, 4] leaves profile distances 9.529 (`u_p`, best coordinate `d_p`) and 8.380 (`u_q`, best coordinate `d_q`). Both sit above the profile budget 3.841. The gain forcing `u_h` is exactly a rescaling of κ by 3.500, which lies inside the interval, and the distance falls to a numerical zero. The readout forcing `u_w` leaves distance 0 on the reduced chart and 5.295 on the full chart that still records W. Nineteen calls are issued. Fifteen are refused. Four are admitted. The SHA-256 of kinetic Θ is `b6133b85bccf3837499979b2ea74e97add32f20cf70e4a71c7fc79d1ed6e36ad` before the calls and the same string after the refusals and the admissions.

No number is taken from either parent deposit's results file. The catalogue is synthetic. An admitted schedule is not a dose and not an efficacy.

This is research only. It is not a medical device, not clinical decision support, not a dose, and not a cure. No document DOI is registered.

See [DISCLAIMER.md](DISCLAIMER.md). The manuscript is [THESIS.md](THESIS.md).

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Manuscript (Chapters 1 to 5, Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata, no document DOI |
| `DISCLAIMER.md` | Research-only boundary |
| `sim/ledger_reduction.py` | Seeded admission ledger, ranks, and leakage distances (seed 20260921; no RNG draws) |
| `sim/candidates.json` | Synthetic catalogue. SHA-256 pinned in the script |
| `sim/results.json` | Numbers cited in Chapter Four |
| `sim/figures/` | Spectra, ledger, distances, steady signatures |

## Reproduce

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/ledger_reduction.py
python3 build_pdf.py
```

NumPy, SciPy and Matplotlib are required for the toy. The PDF step also needs the `markdown` and `weasyprint` packages. Regenerating the script rewrites `sim/results.json` and `sim/figures/`. The candidate-file SHA-256 is pinned inside `sim/ledger_reduction.py`. A byte edit that does not update the pin raises.

## Cite

Ogbonna KE. Admitted known forcings that remain distinguishable after a stiff–sloppy reduction [Internet]. Thesis #28 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-28-admitted-forcings-after-reduction

Machine-readable fields are in `CITATION.cff`. Add a document DOI there only after one exists.

Hub index, for cataloguing only: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

## Licence

Text and sketch code are MIT, with attribution. Computational research only.
