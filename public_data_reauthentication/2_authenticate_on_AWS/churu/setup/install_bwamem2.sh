wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/Source_code_including_submodules.tar.gz
tar zxvf Source_code_including_submodules.tar.gz
cd bwa-mem2-2.2.1
make -j 8
if [ $? -eq 0 ]; then
    ln -srf bwa-mem2* ../../
fi
