# Nanodecoder
c++ tool for identifying barcodes from nanopore reads with a whitelist. 

## troubleshooting
```bash
error while loading shared libraries: libhts.so.2: cannot open shared object file: No such file or directory
```

```bash
cd /usr/lib
sudo wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
sudo tar -vxjf htslib-1.9.tar.bz2
cd htslib-1.9
sudo make
sudo ln -s libhts.so.2 ../
```
