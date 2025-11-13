---
title: Intel Compilers oneAPI
tags:
    - Intel
    - oneAPI
---


## Installation of Intel oneAPI compiler suite (v2022.3.1)

### Installation without sudo access

Download Intel oneAPI basekit and HPCkit for free from

[Intel oneAPI Basekit](https://registrationcenter-download.intel.com/akdlm/irc_nas/18970/l_BaseKit_p_2022.3.1.17310_offline.sh){: .btn .btn-green .mr-4}

[Intel oneAPI HPC kit](https://registrationcenter-download.intel.com/akdlm/irc_nas/18975/l_HPCKit_p_2022.3.1.16997_offline.sh){: .btn .btn-green .mr-4}


```bash
chmod a+x ./l_BaseKit_p_2022.3.1.17310_offline.sh
sh ./l_BaseKit_p_2022.3.1.17310_offline.sh

chmod a+x ./l_HPCKit_p_2022.3.1.16997_offline.sh
sh ./l_HPCKit_p_2022.3.1.16997_offline.sh
```

After installation export the path to your `~/.bashrc` file


### Installation with sudo access

```bash
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
sudo add-apt-repository "deb https://apt.repos.intel.com/oneapi all main"
sudo apt-get update
```

### Install the components

```bash
sudo apt-get install -y intel-oneapi-common-vars
sudo apt-get install -y intel-oneapi-compiler-fortran-2022.3.1
sudo apt-get install -y intel-oneapi-mkl-2022.3.1
sudo apt-get install -y intel-oneapi-mkl-devel-2022.3.1
sudo apt-get install -y intel-oneapi-mpi-2022.3.1
sudo apt-get install -y intel-oneapi-mpi-devel-2022.3.1
```

