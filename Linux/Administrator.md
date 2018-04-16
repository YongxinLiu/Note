<!-- TOC -->

- [1. 用户User](#1-用户user)
    - [1.1.](#11)
- [2. 软件Softwares](#2-软件softwares)
    - [2.1. BioConda](#21-bioconda)
- [3. 数据库Databases](#3-数据库databases)

<!-- /TOC -->

# 1. 用户User

## 1.1. 


# 2. 软件Softwares

## 2.1. BioConda

    #官网https://bioconda.github.io/
    cd ~/Downloads
    wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    
    #确认许可协议，默认安装目录为conda，不添加环境变量
    bash Miniconda3-latest-Linux-x86_64.sh # yes /conda no 
    
    # 手动设置3个常用命令到环境，添加目录会替换为Python3环境
    ln /conda/bin/conda /usr/local/bin/
    ln /conda/bin/activate /usr/local/bin/
    ln /conda/bin/deactivate /usr/local/bin/
    
    #添加bioconda，越靠后优先越高
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/ 
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/ 
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/ 
    conda config --set show_channel_urls yes
    conda config --add channels r # Optional
    
    #显示已有的通道
    conda config --get channels


# 3. 数据库Databases