wget -i Nature2019data/metagenomics_address.txt -P Nature2019data/MGX/
tar -zxvf Nature2019data/MGX/*.tar
cd Nature2019data/MGX_unzip
ls ../MGX/*.tar.bz2 |xargs -n1 tar xvf


wget -i Nature2019data/metatranscriptomics_address_core.txt -P Nature2019data/MTX/
cd Nature2019data/MTX_unzip
ls ../MTX/*.tar.bz2 |xargs -n1 tar xvf
rm *log *_ecs.tsv *pathabundance.tsv
