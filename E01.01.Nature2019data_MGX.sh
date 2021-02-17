wget -i Nature2019data/metagenomics_address.txt -P Nature2019data/MGX/
tar -zxvf Nature2019data/MGX/*.tar
cd Nature2019data/MGX_unzip
ls ../MGX/*.tar.bz2 |xargs -n1 tar xvf

