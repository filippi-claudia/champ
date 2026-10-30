---
title: "CMake"
---

# CMake Installation

!!! note

    CHAMP requires CMake version 3.17 or higher for building.

## Installing CMake

### Using Package Manager

Most modern Linux distributions provide recent versions of CMake through their package managers:

**Ubuntu/Debian:**
```
sudo apt-get install -y cmake
```

**Fedora/RHEL:**
```
sudo dnf install cmake
```

**macOS (Homebrew):**
```
brew install cmake
```

### Installing from Binary Distribution

If your system's package manager provides an older version, you can install the latest CMake from the official releases:

1. Visit the [CMake releases page](https://github.com/Kitware/CMake/releases)
2. Download the appropriate binary for your system (e.g., `cmake-X.Y.Z-linux-x86_64.sh`)
3. Extract and add to your PATH:

```bash
wget https://github.com/Kitware/CMake/releases/download/vX.Y.Z/cmake-X.Y.Z-linux-x86_64.sh
chmod +x cmake-X.Y.Z-linux-x86_64.sh
./cmake-X.Y.Z-linux-x86_64.sh --skip-license --prefix=$HOME/.local
export PATH=$HOME/.local/bin:$PATH
```

Replace `X.Y.Z` with the desired version number (3.17 or higher).

### Building from Source

For the latest features or custom installations:

```bash
wget https://github.com/Kitware/CMake/releases/download/vX.Y.Z/cmake-X.Y.Z.tar.gz
tar -xzvf cmake-X.Y.Z.tar.gz
cd cmake-X.Y.Z
./bootstrap --prefix=$HOME/.local
make -j$(nproc)
make install
export PATH=$HOME/.local/bin:$PATH
```

## Verifying Installation

Check your CMake version:

```
cmake --version
```

Ensure the version is 3.17 or higher.

## HPC Systems

On HPC systems, CMake is often available through the module system:

```bash
module load cmake
```

or

```bash
module load CMake
```

Check available versions with `module avail cmake`.


