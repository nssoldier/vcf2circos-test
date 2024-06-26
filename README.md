# Installation

## local
```
wget https://www.lbgi.fr/~lamouche/vcf2circos/config_vcf2circos_29032023.tar.gz
tar -xzf config_vcf2circos_29032023.tar.gz
mv Static ~/
python setup.py install
python main.py
```
## Docker
```
docker image build -t vcf2circos . --no-cache
docker run -d --restart=no --expose 5000 -p 5000:5000 vcf2circos
```