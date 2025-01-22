git clone https://github.com/lh3/bwa.git
cd bwa
make -j 64
if [ $? -eq 0 ]; then
    ln -srf ./bwa ../../
fi
