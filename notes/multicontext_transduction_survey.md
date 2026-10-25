# Public multi-context single-cell CRISPR data sets: one transduction or separate transductions?

Compiled 2026-10-25.

**Why.** When the input context and the scoring context of a transfer model come from one
transduced and expanded cell population, the two contexts contain descendants of the same
transduced cells (shared clones). Scores measured across such contexts are not a clean test of
transfer. This file records, for each commonly used multi-context single-cell CRISPR data set,
whether the contexts come from ONE transduction or from SEPARATE transductions, with the
supporting methods sentence quoted.

**Method.** Full texts were read from Europe PMC / PMC (`fullTextXML`, NCBI `efetch` fallback)
for every paper with an open-access full text; the transduction, screen-design and
single-cell-collection paragraphs of the methods were read directly. Calls are based on methods
text only; where methods text was not reachable the call is marked UNCLEAR and the sources that
were visible are named. bioRxiv HTML/PDF returns HTTP 403 to automated fetches and was not
bypassed; for the one preprint affected (X-Atlas/Orion) only the author-written data-repository
description was available. Benchmark usage was verified in the benchmark papers / repositories
themselves (see "Benchmarks checked" below), not assumed.

**Call vocabulary**

- **ONE** - all contexts of the relevant axis descend from a single transduced, selected and
  expanded population that was later split, treated, differentiated, transferred or sampled.
  Clone sharing across contexts is expected.
- **SEPARATE** - each context was transduced independently (different cell line, different donor,
  different screen). No clone sharing possible.
- **MIXED** - ONE along one context axis (e.g. treatment, time point) and SEPARATE along another
  (e.g. donor, cell line) in the same data set.
- **UNCLEAR** - methods text not reachable or not explicit.

## Table

| # | Data set / reference (DOI) | Contexts | One or separate transduction | Cells transduced per guide (as stated) | Used in model benchmarks | Quote |
|---|---|---|---|---|---|---|
| 1 | Replogle et al. 2022, *Cell* 185:2559 (10.1016/j.cell.2022.05.013) | K562-genome-wide (day 8), K562-essential (day 6), RPE1-essential (day 7) | **SEPARATE** - three independently titred and transduced experiments; the two K562 data sets are separate spinfections of the same parental CRISPRi line | "coverage of 1000 cells per library element" maintained; median cells/perturbation after filtering 171 (K562-GW), 124 (K562-ess), 79 (RPE1) | scPerturb (Peidli et al. 2024); scPerturbBench (Li et al. 2026); State / Replogle-Nadig (Adduri et al. 2026); PIE (`PIE_replogle_nadig_essential`) | Q1a-c |
| 2 | Nadig et al. 2025, *Nat Genet* (10.1038/s41588-025-02169-3) | Jurkat-essential, HepG2-essential (plus the three Replogle contexts re-analysed) | **SEPARATE** - Jurkat and HepG2 each transduced in their own experiment with library dJR092 | ">1000 cells per library element" each; median cells/perturbation 83 (Jurkat), 45 (HepG2) (Table 1) | State: held-out context `replogle.jurkat` (zero-shot), `replogle.k562` (few-shot); PIE; Zhu et al. 2026 used Jurkat-essential as a cell-line comparison | Q2a-b |
| 3 | Zhu et al. 2026, *Cell* (10.1016/j.cell.2026.08.002); preprint 2025 (10.64898/2025.12.23.696273) | Rest, Stim8hr, Stim48hr x 4 donors; genome-scale CRISPRi in primary CD4+ T cells | **MIXED** - ONE within a donor: one transduced, selected, expanded culture split into the three conditions on day 12. SEPARATE across the four donors | library MOI ~0.2; 30 million naive CD4+ T cells seeded per donor at day 0; per-perturbation coverage shown in Fig. S1C, not stated numerically in text | Not found in the benchmarks checked (data released 2025-12 / 2026-08) | Q3a-b |
| 4 | Frangieh et al. 2021, *Nat Genet* 53:332 (10.1038/s41588-021-00779-1) | control, IFN-gamma, TIL co-culture (Perturb-CITE-seq) | **ONE** - a single ICR-library transduction at MOI 0.1, cultured 14 days, then seeded into the three conditions | "~1,000 cells/guide" at transduction; ">1,800x representation" per condition at seeding | PerturBench (Frangieh21: train on simple conditions, predict co-culture); scPerturb; scPerturbBench | Q4a-b |
| 5 | Jiang et al. 2025, *Nat Cell Biol* (10.1038/s41556-025-01622-z); preprint 2024 (10.1101/2024.01.29.576933) | 6 cell lines (A549, MCF7, HT29, HAP1, BxPC3, K562) x 5 cytokine stimuli = 30 samples | **SEPARATE** - each cell line was infected separately, per pathway library; each sample carries its own stimulus | Not stated in the accessible text (supplementary methods not retrieved) | PerturBench (Jiang24: 219 perturbations, 6 cell lines, 5 cytokine treatments) | Q5 |
| 6 | Shifrut et al. 2018, *Cell* 175:1958 (10.1016/j.cell.2018.10.024) | resting vs stimulated CD8+ T cells, 2 donors (SLICE + CROP-seq) | **MIXED**, partly UNCLEAR - ONE transduced pool per donor carried through Cas9 electroporation, puromycin selection and sorting before restimulation; SEPARATE across donors. The methods describe only restimulation of that pool, not a separately handled unstimulated arm, although Results analyse "resting and stimulated" cells | Not stated for the CROP-seq arm (48 guides) | scPerturb (ShifrutMarson2018, GSE119450) | Q6a-c |
| 7 | Schmidt et al. 2022, *Science* 375:eabj4008 (10.1126/science.abj4008) | stimulated vs nonstimulated primary human T cells, 2 donors (CRISPRa Perturb-seq) | **MIXED** - ONE per donor: one library transduction at MOI 0.3, then both conditions sorted from it; SEPARATE across donors | library MOI 0.3; 80,000 cells loaded per condition | Not found in the benchmarks checked | Q7a-b |
| 8 | Papalexi et al. 2021, *Nat Genet* 53:322 (10.1038/s41588-021-00778-2) | THP-1 stimulated (DAC + IFN-gamma + TGF-beta1) vs unstimulated (ECCITE-seq) | **ONE** per transduction replicate - three independent low-MOI transductions; each pool was then stimulated or left unstimulated. (Unstimulated cells were 10% of the loaded pool and were excluded from downstream analysis) | Not stated; 22,606 cells with a detected guide over 111 guides | scPerturb (PapalexiSatija2021); scPerturbBench | Q8a-b |
| 9 | Tian et al. 2019, *Neuron* 104:239 (10.1016/j.neuron.2019.07.014) | iPSC vs day-7 iPSC-derived neuron (CROP-seq, 58 guides) | **ONE** - one iPSC transduction; after selection and expansion the same pool was either kept as iPSCs or differentiated into neurons | ~20,000 iPSCs and ~20,000 day-7 neurons captured over 58 guides (~340 cells/guide/context) | scPerturb (TianKampmann2019, GSE152988) | Q9 |
| 10 | Tian et al. 2021, *Nat Neurosci* 24:1020 (10.1038/s41593-021-00862-0) | CRISPRi neurons (184 genes) and CRISPRa neurons (100 genes), day-10 neurons only | **SEPARATE** between the CRISPRi and CRISPRa data sets (different effector lines, different libraries); each data set is a single context, so there is no within-data-set context axis | library MOI 0.1-0.2; ~98,000 CRISPRi neurons collected | scPerturb / scPerturbBench (TianKampmann2021_CRISPRi, _CRISPRa) | Q10 |
| 11 | Datlinger et al. 2017, *Nat Methods* 14:297 (10.1038/nmeth.4177) | Jurkat anti-CD3/CD28-stimulated vs unstimulated (naive) (CROP-seq) | **ONE** - a single transduced flask grown 10 days under selection, then serum-starved and split into stimulated and continued-starvation arms | Not stated per guide | scPerturb (DatlingerBock2017); scPerturbBench Task 1 | Q11a-b |
| 12 | Dixit et al. 2016, *Cell* 167:1853 (10.1016/j.cell.2016.11.038) | BMDC 0 h vs 3 h post-LPS; K562 TF pool day 7 vs day 13; K562 low vs high MOI | **ONE** on both the LPS and the time-point axes - one transduced population per pool, harvested at two times | Table S1 (per-pool); K562 low-MOI culture 10% BFP+, high-MOI 50% BFP+ | scPerturb (DixitRegev2016); widely used in GEARS/CPA-era work | Q12a-b |
| 13 | Zhou et al. 2023, *Nature* 624:154 (10.1038/s41586-023-06733-x) | in vitro input aliquot vs intratumoral (day 7) CTLs; 60 recipient mice | **ONE** - one transduced, Ametrine-sorted OT-I population; an aliquot kept as input and 4x10^6 cells transferred into the mice. The single-cell screen itself reads one tissue and one time point | 4x10^6 transduced cells transferred over 180 targets (dual-guide library); 1x10^6 kept as input | Not found in the benchmarks checked | Q13a-b |
| 14 | Knudsen et al. 2025, *Nature* (10.1038/s41586-025-09489-8) | bulk screen: end-of-production (in vitro), day 7 and day 21 bone marrow. Perturb-seq: day-21 bone marrow only | **ONE** per donor product across the three bulk contexts (one sequential CAR + library transduction, then split in vitro / into mice); SEPARATE across healthy donors. The Perturb-seq arm is single-context | CAR and sgRNA library each at MOI 5; per-guide cell numbers not stated | Not found in the benchmarks checked | Q14a-c |
| 15 | Fagerberg et al. 2025, *Science* 389:eadn2337 (10.1126/science.adn2337) | in vivo perturbSEQ: spleen, day 8 of acute LCMV - single context (one transduced P14 pool across recipient mice) | **ONE** across recipient mice; no second perturbed context in the single-cell data | low MOI stated, no number; 59,815 cells over ~39 targets x 3 guides | Not found in the benchmarks checked | Q15a-b |
| 16 | McFaline-Figueroa et al. 2024, *Cell Genomics* 4:100487 (10.1016/j.xgen.2023.100487) | 3 glioblastoma lines (A172, T98G, U87MG) x drugs/doses (lapatinib, nintedanib, trametinib, zstk474, temozolomide, 6-thioguanine) | **MIXED** - ONE within a cell line: a transduced, puromycin-selected pool was seeded into the drug/dose wells; SEPARATE across cell lines, and the kinome screen spans "two independent transductions" | library MOI ~0.1; 3x10^6 cells transduced per screen | PerturBench (McFalineFigueroa23: 3 cell lines, 525 genetic perturbations, 5 chemicals) | Q16a-c |
| 17 | X-Atlas/Orion, Huang et al. 2025, bioRxiv (10.1101/2025.06.11.659105); data 10.25452/figshare.plus.29190726 | HCT116 and HEK293T, genome-wide CRISPRi | **SEPARATE** across the two lines by construction (two genome-wide experiments in two different lines); **UNCLEAR** otherwise - the methods text could not be read (bioRxiv returns HTTP 403 to automated fetches; no PMC/Europe PMC full text). What was visible: the authors' figshare record, the Hugging Face copy, and the dataset's own guide library table | Not visible in the accessible sources | Q17 |
| 18 | Virtual Cell Challenge 2025 H1 data; Roohani et al. 2025, *Cell* (10.1016/j.cell.2025.06.008) | H1 human embryonic stem cells, CRISPRi - single training context (held-out perturbations, not held-out contexts) | **UNCLEAR** - no methods text reachable: the Cell article page returns HTTP 403 and virtualcellchallenge.org serves a JavaScript shell with no dataset text. The published abstract describes "purpose-built datasets" without protocol detail, so no transduction call is made here. Since the released data are one context, the within-/across-context distinction does not arise for it | Not stated in any accessible source | State (VCC training recipe), PIE (`PIE_arc_vcc_25`); VCC 2026 (10.1016/j.cell.2026.08.004) moves to unseen cell lines | Q18 |

Not examined here (single-context K562 data sets with no multi-context axis, listed only so their
absence is not read as an omission): Adamson et al. 2016, Norman et al. 2019.

## Quotes

Each quote is from the paper's own methods (or, where marked, results) text; section names are the
paper's own headings.

**Q1a** Replogle et al. 2022, STAR Methods, "Perturb-seq experiments" (K562 day 8 genome-scale):
"CRISPRi K562 cells were transduced via spinfection (1000g) with polybrene (8 ug/ml) with the
target of obtaining an infection rate of ~10%."

**Q1b** Same section, K562 day 6 essential-scale: "To perform our K562 day 6 essential-scale
Perturb-seq experiment, library lentivirus was packaged into lentivirus in 293T cells and
empirically measured in K562 cells to obtain viral titers." (The transduction sentence quoted in
Q1a appears once in each of the two K562 paragraphs - one spinfection per experiment.)

**Q1c** Same section, RPE1 day 7: "CRISPRi RPE1 cells expressing ZIM3 KRAB-dCas9-P2A-BFP were
transduced via replating into virus-laden media with polybrene (8 ug/ml) with the target of
obtaining an infection rate of ~10%." Coverage, same section: "Cells were maintained at a
viability of >90%, a coverage of 1000 cells per library element". Per-perturbation medians are in
the Figure S1F-H legends.

**Q2a** Nadig et al. 2025, Methods, "Perturb-seq experiments, library preparation, and
sequencing": "For the Jurkat Perturb-seq experiment, Jurkat cells expressing Zim3-dCas9-P2A-mCherry
were transduced with dJR092 library lentivirus by spinfection (1000g) with polybrene".

**Q2b** Same section: "For the HepG2 Perturb-seq experiment, HepG2 cells expressing dCas9-BFP-KRAB
were transduced with dJR092 library lentiviral supernatant with polybrene (8 ug/mL, Sigma-Alrich)
with a targeted low infection rate of ~10%."

**Q3a** Zhu et al. 2026, STAR Methods, "Genome-scale perturb-seq": "At day 12 morning, cells were
centrifuged and resuspended at 1 x 106 cells/ml with cXVIVO supplemented with 200 IU/mL IL-2 and
split into three populations."

**Q3b** Same section: "In the afternoon, cells were transduced with 0.1% v/v (multiplicity of
infection ~ 0.2) genome-scale gRNA library lentivirus." Donor structure, Results, "A scalable
platform for perturb-seq in primary CD4+ T cells": "Primary naive CD4+ T cells isolated from four
human donors were stimulated, transduced with CRISPRi machinery and the genome-scale library, and
expanded ex vivo".

**Q4a** Frangieh et al. 2021, Results: "cultured the cells for 14 days, and then either
co-cultured them with a range of TIL doses (1:1, 2:1 or 4:1 ratios), treated them with IFNgamma
(no co-culture), or maintained them in culture media alone (control)."

**Q4b** Methods, "Perturb-CITE-Seq co-culture experiment": "After 14 days of editing and
selection, 1.25x105 Cas9-expressing, ICR-library perturbed melanoma target cells were seeded per
well of a 12-well plate (One Plate with 1.5x106 cells total per condition, >1,800x
representation)." Guide coverage, Results: transduction "at multiplicity of infection (MOI) 0.1 to
achieve a high rate of single-guide transductions with ~1,000 cells/guide".

**Q5** Jiang et al. 2025, Results, "Scalable and flexible Perturb-seq across cell lines and
conditions": "For each pathway, we separately infected all six cell lines with the pathway-specific
sgRNA library, and then stimulated infected cells with the corresponding cytokine for 24 hours to
activate signaling".

**Q6a** Shifrut et al. 2018, STAR Methods, "SLICE adapted to CROP-Seq": "Lentivirus was produced
from this pooled plasmid library and used to transduce CD8+ T cells from two healthy donors, as
above."

**Q6b** Same section, after Cas9 electroporation, puromycin selection and live-cell sorting: "Two
days post sorting, cells were restimulated as above. 36 hours post restimulation, cells were
collected, counted, and prepared for Illumina sequencing".

**Q6c** Results, "SLICE Paired with Single Cell RNA-Seq ...": "we analyzed the transcriptional
states of more than 15,000 resting and stimulated single cells where we could identify an sgRNA
barcode." (The methods text quoted in Q6b does not describe how the resting cells were held back,
hence the partial UNCLEAR.)

**Q7a** Schmidt et al. 2022, Materials and Methods, "Perturb-seq sample preparation and
sequencing": "Bulk CD3+ primary human T cells from two donors were transduced and cultured as
described in the 'Genome-wide CRISPRa and CRISPRi screens' section, except library transduction was
completed at lower MOI of 0.3."

**Q7b** Same section: "Twenty-four hours later, cells from both the stimulated and nonstimulated
condition were sorted for mCherry+ (marking dCas9-VP64)."

**Q8a** Papalexi et al. 2021, Results, "ECCITE-seq validates putative immune checkpoint
regulators": "Our total dataset represents three independent transductions at low multiplicity of
infection, aiming to maximize the proportion of cells infected with a single gRNA."

**Q8b** Methods, "ECCITE-seq experimental setup": "After five days of doxycycline treatment, cells
were stimulated with DAC, IFNgamma and TGFbeta1 for an additional three days or left unstimulated
prior to running the 10x Genomics experiment".

**Q9** Tian et al. 2019, STAR Methods, "CROP-Seq": "CRISPRi-i3 N iPSCs were infected with a pool
of selected sgRNAs in the CROP-Seq vector pMK1334 at a low multiplicity of infection ... After
puromycin selection and expansion, cells were either passaged as iPSCs or differentiated into
neurons."

**Q10** Tian et al. 2021, Methods (CROP-seq): "Day 0 CRISPRi and CRISPRa neurons were infected by
the corresponding CROP-seq sgRNA library at a MOI of 0.1-0.2, followed by puromycin selection at 4
ug/ml for 3 days and recovery."

**Q11a** Datlinger et al. 2017, Online Methods, "Lentiviral transduction with gRNA libraries or
single gRNAs": "A flask with about 30% surviving cells was chosen and grown for 10 d to allow for
efficient genome editing while renewing the selective medium".

**Q11b** Online Methods, "Anti-CD3/CD28 stimulation of the T cell receptor pathway in Jurkat
cells": "They were stimulated with 1 ug/ml anti-CD3 ... and 1 ug/ml anti-CD28 antibody ... for 4 h,
while the unstimulated (naive) control was subjected to continued starvation. Samples collected for
both conditions were subjected to CROP-seq."

**Q12a** Dixit et al. 2016, Results (BMDC): "After another 7 days, we stimulated the cells with
LPS, and collected cells for scRNA-seq at 0 and 3h (32,624 and 37,369 cells, respectively, Table
S1)."

**Q12b** Methods and Resources, "K562 cell cultures": "Cells were harvested for library preparation
7 days post transduction for most experiments and 13 days post transduction for the second time
point of the TF pool experiment."

**Q13a** Zhou et al. 2023, Methods, "In vivo screening": "Transduced cells were sorted based on the
expression of Ametrine, and an aliquot of 1 x 106 transduced OT-I cells was saved as input."

**Q13b** Same section: "Transduced OT-I cells (4 x 106) were then transferred intravenously to
B16-OVA tumour-bearing C57BL/6 mice at day 12 after tumour inoculation. A total of 60 recipient
mice was used in 2 experiments combined."

**Q14a** Knudsen et al. 2025, Methods, "Lentivirus and CAR-T cell production": "For production of
Mario-CAR-T cells, activated cells were transduced sequentially first with the BCMA CAR and then
with Mario sgRNA library lentivirus (both MOI 5)."

**Q14b** Methods, "In vivo loss of function CRISPR screen": "Mice were euthanized at either day 7
or day 21, and femur, tibia, and spine were collected for isolation of marrow infiltrating T
cells."

**Q14c** Results, "In vivo Perturb-seq of BCMA CAR-T cells": "Following transduction with the
sgRNA library and CAR vectors, engineered T cells were enriched by negative magnetic bead selection
and transferred into mice previously engrafted with MM.1S myeloma". The single-cell readout is one
context: "After 21 days, mice were euthanized and the modified CAR-T cells were [isolated] from the
bone marrow".

**Q15a** Fagerberg et al. 2025, Results, "Performing perturbSEQ on select TFs and EMs in acute
LCMV infection": "dsRED+ Cas9+ P14 CD8 T cells were transduced at a low multiplicity of infection
(MOI), positively selected on hCD2, and transferred into day 1 acute LCMV infected congenic
recipients".

**Q15b** Same section: "At day 8 post infection, spleens were harvested and P14s and endogenous
H2Db-gp33 specific CD8 T cells were sorted".

**Q16a** McFaline-Figueroa et al. 2024, STAR Methods, "Generation of CROP-seq-OPTI gRNA libraries":
"For our screens, 3 x 106 cells in 10 cm tissue culture dishes were transduced with lentiviral
libraries at an approximate MOI of 0.1 to ensure single integrations."

**Q16b** STAR Methods, "Chemical exposure of genetically perturbed cell pools": "Genetically
perturbed pools of glioblastoma cells were seeded in 96-well plates at 2.5 x104 cells per well ...
and 1 uL added of the appropriate drug and dose to wells of seeded cells".

**Q16c** Results: "Our screen contained 14,121 gene-by-environment combinations across two
independent transductions for 28,242 unique conditions across 1,052,205 single-cell
transcriptomes".

**Q17** X-Atlas/Orion: no methods text was reachable. The authors' data record states: "This
dataset (X-Atlas/Orion) contains processed data from two genome-wide Perturb-seq experiments in
HCT116 and HEK293T cell lines" (figshare+ item 29190726, dataset description). The guide library
table in the same record states that guide sequences are from Replogle et al., eLife (2022). No
transduction, MOI or per-guide coverage statement was visible.

**Q18** Virtual Cell Challenge 2025: no methods text was reachable. The published abstract states
only that the challenge "will provide an evaluation framework, purpose-built datasets, and a venue
for accelerating model development" (Roohani et al. 2025, Cell, abstract). The 2026 edition's
abstract states its new data comprise "unseen cell lines" (Cell 2026, 10.1016/j.cell.2026.08.004,
abstract), which is consistent with the 2025 release being a single context, but gives no
transduction detail.

## Benchmarks checked

- **scPerturb** - Peidli et al. 2024, *Nat Methods* 21:531 (10.1038/s41592-023-02144-y). Harmonized
  resource; the data-availability list includes DatlingerBock2017, DixitRegev2016, FrangiehIzar2021,
  PapalexiSatija2021, ReplogleWeissman2022, ShifrutMarson2018 (GSE119450), TianKampmann2019
  (GSE152988), TianKampmann2021 (GSE124703).
- **scPerturbBench** - Li et al. 2026, "A systematic comparison of single-cell perturbation
  response prediction models", *Sci Adv* (10.1126/sciadv.aed3414; code
  github.com/TianGzlab/scPerturbBench). Its dataset table includes Datlinger, Adamson, Frangieh,
  Papalexi, Replogle (RPE1), Tian 2021 CRISPRi and CRISPRa, Norman, Srivatsan, and a cell-state
  transfer task.
- **PerturBench** - Wu et al., arXiv:2408.10609 (NeurIPS 2024/2025 Datasets & Benchmarks). Table 1
  datasets: Srivatsan20, Frangieh21 (248 genetic perturbations; train on the simpler conditions,
  predict co-culture), Jiang24 (219 perturbations, 6 cell lines, 5 cytokine treatments),
  McFalineFigueroa23 (3 cell lines, 525 genetic perturbations, 5 chemicals), Norman19, OP3.
- **State** - Adduri et al. 2026, *Cell* (10.1016/j.cell.2026.07.052); preprint
  10.1101/2025.06.26.661135; code github.com/ArcInstitute/state. Its example configurations train on
  the Replogle-Nadig dataset with held-out contexts `replogle.jurkat` (zero-shot) and
  `replogle.k562` (few-shot), plus Tahoe-100M and the Virtual Cell Challenge recipe.
- **PIE** - Verma et al. 2026, preprint (10.64898/2026.10.02.756297). Public dataset copies
  `PIE_replogle_nadig_essential`, `PIE_tahoe100m`, `PIE_x_atlas_orion`, `PIE_arc_vcc_25`
  (huggingface.co/arcinstitute).
- **Virtual Cell Challenge** - Roohani et al. 2025, *Cell* (10.1016/j.cell.2025.06.008); Arc Virtual
  Cell Initiative Team 2026, *Cell* (10.1016/j.cell.2026.08.004).

## What this means for the preprint's claim

Transfer evaluations that cross contexts **within** one of the ONE/MIXED data sets above
(Frangieh control -> IFN-gamma -> co-culture; Papalexi unstimulated -> stimulated; Datlinger naive
-> stimulated; Dixit 0 h -> 3 h and day 7 -> day 13; Tian 2019 iPSC -> neuron; Shifrut resting ->
stimulated within a donor; Zhu Rest -> Stim8hr -> Stim48hr within a donor; Schmidt nonstimulated ->
stimulated within a donor; McFaline-Figueroa vehicle -> drug within a cell line) score a model on
descendants of the same transduced cells that produced its input. Transfers that cross cell lines
or donors (Replogle K562 -> RPE1; Replogle -> Nadig Jurkat/HepG2, i.e. the State zero-shot task;
Jiang across its six lines; X-Atlas/Orion HCT116 vs HEK293T; any across-donor split in Zhu,
Schmidt, Shifrut or Knudsen) cannot share clones.
