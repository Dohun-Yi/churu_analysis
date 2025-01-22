# x64 architecture
architecture=$(uname -m)

cd $(dirname $0)
if [ "$architecture" == "x86_64" ]; then
    """ x86_64 """
    wget https://github.com/refresh-bio/KMC/releases/download/v3.2.2/KMC3.2.2.linux.x64.tar.gz
    mkdir -p ./KMC/
    tar zxvf KMC3.2.2.linux.x64.tar.gz -C ./KMC/
    ln -srf ./KMC/bin/* ../
else
    """ arm64 """
    wget https://github.com/refresh-bio/KMC/releases/download/v3.2.2/KMC3.2.2.linux.arm64.tar.gz
    mkdir -p ./KMC/
    tar zxvf KMC3.2.2.linux.arm64.tar.gz -C ./KMC/
    ln -srf ./KMC/bin/* ../
fi

## install from source
#git clone --recurse-submodules https://github.com/refresh-bio/KMC.git
#cd KMC/
#./build_release.py
