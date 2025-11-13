---
layout: default
title: CMake
parent: Dependencies
nav_order: 1
authors:
    - Ravindra Shinde
tags:
    - CHAMP
    - CMake
---


## Install or load [cmake](https://cmake.org/)

```bash
sudo apt-get install -y cmake
```

or

```bash
wget https://github.com/Kitware/CMake/releases/download/v3.25.1/cmake-3.25.1-linux-x86_64.sh
tar -xzvf cmake-3.25.1-linux-x86_64.sh
export PATH=cmake-3.25.1-linux-x86_64/bin:$PATH
```

{: .warning }
The version of the cmake must be greater than 3.17.


