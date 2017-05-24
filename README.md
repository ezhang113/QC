# QC summary parser

### Items Parsed for RNA-SEQ

- *.NRF.metrics <i>unused</i> : NRF
- *.fqTr.metrics : first cutadapt trim
- *.REPELEMENT.metrics : repetitive mapped elements file <i>coming soon</i>
- *.fqTr\*U-SoMa.metrics : STAR map log (2nd mapping)

### Items parsed for eCLIP
- *.fqTrTr.metrics : second cutadapt trim
- *.fqTrTrU-SoMaSoCo.metrics : UMI-duplicate removed metrics
- *.fqTrTrU-SoMaSoCoSoV2Cl.bed.log <i>unused</i>
- *.fqTrTrU-SoMaSoCoSoMeV2Cl.bed : clipper peaks file
- *.fqTrTrU-SoMaSoCoSoMeV2ClNoCoFc1Pv1.204_01.metrics <i> would be useful to have at some point </i>

# TODO:
1. add repetitive element stats (ie. number/percentage of reads mapping to rRNA)
   - use: data/old/PKM2_IN_R1_M_S17_L001_R1_001.unassigned.adapterTrim.round2.rmRep.metrics
2. look into adding the input-normalized peak counts to the metrics file
3. the "204_01" is being added to input normalized peaks, is this an artifact/hardcoded somewhere?
