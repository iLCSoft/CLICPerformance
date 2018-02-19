# Flavour tagging

## Description

Folder containing steering files to perform the vertex reconstruction, jet clustering and produce flavour-tagging weights from samples with a given number of jets.
Additionally it contains a python script to submit an example job and a script to evaluate the flavour tagging efficiencies.
The flavour tagging analysis is based on the [LCFIPlus](https://github.com/lcfiplus/LCFIPlus) framework within iLCSoft. For more information, please go to https://doi.org/10.1016/j.nima.2015.11.054.

## Content and usage

1. [vtx_makentp.xml](https://github.com/nachogargar/CLICPerformance/blob/add_flavtag_files/Flavour_tagging/vtx_makentp.xml): Steering file that performs the vertex reconstruction, the jet clustering and produce ntuples with the required information for the flavour-tagging training (impact parameters, vertex and jet properties).
- **Input**: output LCIO file from [clicReconstruction.xml](https://github.com/iLCSoft/CLICPerformance/blob/master/examples/clicReconstruction.xml)
- **Output**: 1. LCIO file with MCParticle, PandoraPFOs and Primary/Secondary Vertex information; 2. Ntuple ready for the training. *Repeat this step for each jet flavour (b, c and light flavour)*

2. [training.xml](https://github.com/nachogargar/CLICPerformance/blob/add_flavtag_files/Flavour_tagging/training.xml): Steering file that reads the jet-flavour information ntuples previously created (b, c and light flavour jets) and carry out the training to produce the weight files for each flavour.
- **Input**: 3 root files: Output_Ntuple_bb.root, Output_Ntuple_cc.root, Output_Ntuple_qq.root
- **Output**: folder with the flavour-tagging weights

3. [Job_example.py](https://github.com/nachogargar/CLICPerformance/blob/add_flavtag_files/Flavour_tagging/Job_example.py): Python script to submit a job to ILCDirac

4. [Flavtag_performance.C](https://github.com/nachogargar/CLICPerformance/blob/add_flavtag_files/Flavour_tagging/Flavtag_performance.C): Script that reads the flavour-tagging weights and produce performance plots. This particular example reads two folders (with and without gg->hadrons background) and plots the misidentification efficiency as a function of the b(c)-effciency.
- Within the ROOT environment type the command: **.x Flavtag_performance.C(false,"B")** -> if *false*(*true*) is used the test(train) tree is used. Type *"C"* to plot charm-efficiency.
- **Input**: folder with the flavour-tagging weights
- **Output**: flavour-tagging performance plots


## Some tips

- Path to generated e+e- -> dijet samples: 
> /ilc/user/p/proloff/stdhep/2f_samples_23_04_2013
- When gg->hadrons is overlaid, change the PFO collection to the desired (Loose, Selected, Tight), for instance:
```xml
<parameter name="PFOCollection" type="string" value="SelectedPandoraPFOs" />
```
and uncomment the FasJetProcessor (optimal jet clustering algorithm parameters may vary for different energies)
```xml
<processor name="MyFastJetProcessor"/>
```
- In the steering file [training.xml](https://github.com/nachogargar/CLICPerformance/blob/add_flavtag_files/Flavour_tagging/training.xml) use:
```xml
nvtx==1&&nvtxall==1 for running Marlin locally
nvtx==1&amp;&amp;nvtxall==1 for GRID jobs
```
- For flavour-tagging studies set the parameter *BeamspotConstraint* to 1, for vertex resolution studies use 0:
```xml
<parameter name="PrimaryVertexFinder.BeamspotConstraint" type="bool">1 </parameter>
```
