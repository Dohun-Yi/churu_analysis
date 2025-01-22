git clone https://github.com/alexdobin/STAR.git
cd ./STAR/source/
make STAR CXXFLAGS_SIMD="-march=native -mtune=neoverse-n1 -O3" -j 8
if [ $? -eq 0 ]; then
    ln -srf STAR ../../../
fi
