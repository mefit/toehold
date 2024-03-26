#!/usr/bin/env python3
import json
import numpy as np
import nupack
import sys
# config.parallelism = True

# Specify all strands
if len(sys.argv) < 2:
    raise RuntimeError('No enough inputs specified.')
else:
    inputs = json.loads(sys.argv[1])
    genome = nupack.Strand(inputs[0], name='genome')
    trigger = nupack.Strand(inputs[1], name='trigger')
    switch = nupack.Strand(inputs[2], name='switch')

# Specify the universal model (37°C, 1.0 M Na+, and 0 M Mg2+)
model = nupack.Model(material='rna', ensemble='stacking', celsius=37, sodium=1.0, magnesium=0.0)

# Analyze complexes
genome_alone = nupack.complex_analysis(complexes=[nupack.Complex([genome])], model=model, compute=['pairs', 'mfe'])['(genome)']
switch_alone = nupack.complex_analysis(complexes=[nupack.Complex([switch])], model=model, compute=['pairs', 'mfe'])['(switch)']
target_complex = nupack.complex_analysis(complexes=[nupack.Complex([genome, switch])], model=model, compute=['mfe'])['(genome+switch)']
genome_alone_ptab = genome_alone.pairs.to_array()
switch_alone_ptab = switch_alone.pairs.to_array()

'''favourability (the higher, the better)

A modified version based on (Green et al., 2014) reconsidering the favourable
direction (highest score is the most favourable) and the calculation of the
normalized ensemble defect of the switch (focus on the stem-loop region).
'''

# Local single-strandness of the genome at the sensor binding site
genome_lss = np.sum(
    np.diagonal(genome_alone_ptab)[
        str(genome).find(str(trigger)) + np.arange(len(str(trigger)))
    ]
) / trigger.nt()

# Local single-strandness of the sensor in the toehold switch
switch_lss = np.sum(
    np.diagonal(switch_alone_ptab)[
        str(switch).find(str(~trigger)) + np.arange(len(str(trigger)))
    ]
) / trigger.nt()

# Normalized ensemble defect of the stem-loop in the toehold switch
stemloop_pairlist = nupack.PairList('(11.3(5.12)5.3)11').pairlist()
stemloop_padding = str(switch).find(str(~trigger)) + 25
stemloop_ned = np.sum(
    [switch_alone_ptab[
        stemloop_padding + index, stemloop_padding + stemloop_pairlist[index]
    ] for index in range(stemloop_pairlist.size)]
) / stemloop_pairlist.size

favourability = 5 * genome_lss + 4 * switch_lss + 3 * stemloop_ned

'''delta_mfe (the higher, the better)

Delta MFE between the switch/trigger pair and the switch-trigger complex
'''

mfe_genome = genome_alone.mfe[0].energy
mfe_switch = switch_alone.mfe[0].energy
mfe_complex = target_complex.mfe[0].energy

delta_mfe = mfe_genome + mfe_switch - mfe_complex

# Output
print(json.dumps(dict(
    favourability = favourability,
    delta_mfe = delta_mfe,
    genome_lss = genome_lss,
    switch_lss = switch_lss,
    stemloop_ned = stemloop_ned,
    genome_mfe = mfe_genome,
    genome_mfe_structure = str(genome_alone.mfe[0].structure),
    switch_mfe = mfe_switch,
    switch_mfe_structure = str(switch_alone.mfe[0].structure),
    complex_mfe = mfe_complex,
    complex_mfe_structure = str(target_complex.mfe[0].structure)
)))
