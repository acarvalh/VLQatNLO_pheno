# VLQatNLO_pheno

You need to install https://github.com/delphes <br />
(+ root + matplotlib)

Example of analysis running: <br />

source  /afs/cern.ch/sw/lcg/external/gcc/4.9.3/x86_64-slc6/setup.sh <br />
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.00/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh <br />
(or: load your local enviroment for root and python) <br />

cd $my_repo/VLQatNLO_pheno <br />
./do_links.sh # check the path relative to delphes == run only once <br />
./analysis_gen.py ref_mass sample.root <br />
<br />
'ref_mass' is the mass to consider for doing a specific cut <br />
==> the path for the sample.root is taken relatively to the EOS path from the samples in cernbox == if you are doing this locally take care of the path defined to sample.root on analysis_gen.py <br />
==> If instead of sample.root you add a name without '.root' it will process the QCD samples (taking into acount the relative cross sections)<br />

<br />

# Example of processing a hepmc sample

cd $my_repo/Delphes <br />
source  /afs/cern.ch/sw/lcg/external/gcc/4.9.3/x86_64-slc6/setup.sh <br />
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.00/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh <br />
(or: load your local enviroment for root and python) <br />

cp $my_repo/VLQatNLO_pheno/gen_card.tcl cards <br />
gunzip sample.hepmc.gz  <br />
./DelphesHepMC cards/gen_card.tcl sample.root  sample.hepmc  <br />
gzip sample.hepmc  <br />

There are examples of scripts to ntuplize signal and BKG (process_*)
