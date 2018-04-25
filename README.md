# VLQatNLO_pheno

You need to install https://github.com/delphes \n
(+ root + matplotlib)

Example of analysis running: \n

source  /afs/cern.ch/sw/lcg/external/gcc/4.9.3/x86_64-slc6/setup.sh \n
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.00/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh \n
(or: load your local enviroment for root and python) \n

cd $my_repo/VLQatNLO_pheno \n
./do_links.sh # check the path relative to delphes == run only once \n

\n
########
Example of processing a hepmc sample \n


cd $my_repo/Delphes \n
source  /afs/cern.ch/sw/lcg/external/gcc/4.9.3/x86_64-slc6/setup.sh \n
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.00/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh \n
(or: load your local enviroment for root and python) \n

cp $my_repo/VLQatNLO_pheno/gen_card.tcl cards \n
gunzip sample.hepmc.gz  \n
./DelphesHepMC cards/gen_card.tcl sample.root  sample.hepmc  \n
gzip sample.hepmc  \n
