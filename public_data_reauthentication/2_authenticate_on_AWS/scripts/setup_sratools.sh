git clone https://github.com/ncbi/ncbi-vdb.git
git clone https://github.com/ncbi/sra-tools.git

sudo yum -y install cmake perl-core perl-XML-LibXML libstdc++-static
cd ncbi-vdb/; ./configure; make -j 64; make install; cd -
cd sra-tools/; ./configure; make -j 64; make install; cd -

echo "export PATH=$PATH:/usr/local/ncbi/sra-tools/bin/" >> ~/.bashrc
exec $SHELL -l
