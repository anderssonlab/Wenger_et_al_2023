# Wenger_et_al_2023
Code for analyzing sequencing data used in the Wenger et al. 2023 study.

Pre-processing, mapping and QC of ChIP/SCAR datasets was performed with the ENCODE's official ChIP-seq pipeline:

https://github.com/ENCODE-DCC/chip-seq-pipeline2

Computation of SCAR-seq Partition files was performed with a modified version
of https://github.com/anderssonlab/Replication_SCARseq/blob/master/SCARseq_commands.sh
which includes bias correction from stranded inputs and support for paired-end
data. See SCAR_processing.sh for usage. 

Steps for processing of external datasets can be found in Process_external.R







