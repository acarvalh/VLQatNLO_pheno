# VLQatNLO_pheno

You need to install https://github.com/delphes <br />
(+ root + matplotlib)

Example of analysis running: <br />

source  /afs/cern.ch/sw/lcg/external/gcc/4.9.3/x86_64-slc6/setup.sh <br />
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.00/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh <br />
(or: load your local enviroment for root and python) <br />

cd $my_repo/VLQatNLO_pheno <br />
./do_links.sh # check the path relative to delphes == run only once <br />

<br />
######## <br />
Example of processing a hepmc sample <br />


cd $my_repo/Delphes <br />
source  /afs/cern.ch/sw/lcg/external/gcc/4.9.3/x86_64-slc6/setup.sh <br />
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.00/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh <br />
(or: load your local enviroment for root and python) <br />

cp $my_repo/VLQatNLO_pheno/gen_card.tcl cards <br />
gunzip sample.hepmc.gz  <br />
./DelphesHepMC cards/gen_card.tcl sample.root  sample.hepmc  <br />
gzip sample.hepmc  <br />
