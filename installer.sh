#!/bin/sh

homeDIR="$( pwd )"

echo "Installation will take place in $homeDIR"

echo "[Checking system dependencies]"
PKG_OK=$(dpkg-query -W -f='${Status}' autoconf 2>/dev/null | grep -c "ok installed")
if test $PKG_OK = "0" ; then
  echo "autoconf not found. Install it with sudo apt-get install autoconf."
  exit
fi
PKG_OK=$(dpkg-query -W -f='${Status}' libtool 2>/dev/null | grep -c "ok installed")
if test $PKG_OK = "0" ; then
  echo "libtool not found. Install it with sudo apt-get install libtool."
  exit
fi
PKG_OK=$(dpkg-query -W -f='${Status}' gzip 2>/dev/null | grep -c "ok installed")
if test $PKG_OK = "0" ; then
  echo "gzip not found. Install it with sudo apt-get install gzip."
  exit
fi
PKG_OK=$(dpkg-query -W -f='${Status}' bzr 2>/dev/null | grep -c "ok installed")
if test $PKG_OK = "0" ; then
  echo "bzr not found. Install it with sudo apt-get install bzr."
  exit
fi

cd $homeDIR


#madgraph="MG5_aMC_v3.6.2.tar.gz"
#URL=https://launchpad.net/mg5amcnlo/3.0/3.6.x/+download/$madgraph
# Use github version with fix for printing intermediate particles
madgraph="mg5amcnlo-3.6.3.zip"
echo -n "Install MadGraph (y/n)? "
read answer
if echo "$answer" | grep -iq "^y" ;then
	#mkdir MG5;
	#echo "[installer] getting MadGraph5"; wget $URL 2>/dev/null || curl -O $URL; tar -zxf $madgraph -C MG5 --strip-components 1;
	cp $madgraph $madgraph"_bak"
	unzip $madgraph;
	mv mg5amcnlo-3.6.3 MG5;
	mv $madgraph"_bak" $madgraph
	cd $homeDIR
	cd ./MG5/bin;
	echo "[installer] installing HepMC, LHAPDF6 and Pythia8 under MadGraph5"
        echo "install hepmc\ninstall lhapdf6\ninstall pythia8\nexit\n" > mad_install.txt;
	./mg5_aMC -f mad_install.txt
	cd $homeDIR
	sed  "s|homeDIR|$homeDIR|g" mg5_configuration.txt > ./MG5/input/mg5_configuration.txt;
fi

echo -n "Install Delphes (y/n)? "
read answer
if echo "$answer" | grep -iq "^y" ;then
# Check if pythia8 has been installed under MadGraph5
  pythiaDir=$homeDIR/MG5/HEPTools/pythia8 
  if [ ! -d "$pythiaDir" ]; then
    echo "Delphes should be installed after hepmc, lhapdf6 and pythia8 were installed in MadGraph."
    exit
  fi
  echo "[installer] Installing DelphesLLP";
  
  tar -zxf DelphesLLP.tar.gz;
  cd DelphesLLP;
  export PYTHIA8=$pythiaDir;
  make HAS_PYTHIA8=true;
  rm -rf .git
  cd $homeDIR;
#  rm DelphesLLP.tar.gz;
fi

echo "\n[installer] For running Delphes the following env variables should be set:\n\n export LD_LIBRARY_PATH=$homeDIR/MG5/HEPTools/pythia8/lib"
echo "\nand for reading Delphes produced ROOT files:\n\n export ROOT_INCLUDE_PATH=$homeDIR/DelphesLLP/external\n"
echo "\n\n or run source setenv.sh\n\n"
cd $currentDIR
