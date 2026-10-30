---
title: Installing with Docker
tags:
    - installation
    - docker
    - containers
---

# Installing CHAMP with Docker

CHAMP is available as pre-built container images on [Docker Hub](https://hub.docker.com/r/neelravi/champ), providing a quick and easy way to run CHAMP without compiling from source. Docker containers include all necessary dependencies and are ready to use immediately.

## What is Docker?

Docker is a containerization platform that packages applications and their dependencies into portable containers. Using Docker for CHAMP provides:

- **No compilation required** - Pre-built binaries ready to use
- **Consistent environment** - Same setup across different systems
- **Isolated dependencies** - No conflicts with system libraries
- **Easy updates** - Pull new versions with a single command
- **Multiple configurations** - Different compiler and library combinations available

## Prerequisites

### Install Docker

**Ubuntu/Debian:**
```
sudo apt-get update
sudo apt-get install docker.io
sudo systemctl start docker
sudo systemctl enable docker
```

**Fedora/RHEL:**
```
sudo dnf install docker
sudo systemctl start docker
sudo systemctl enable docker
```

**macOS:**

Download and install [Docker Desktop for Mac](https://www.docker.com/products/docker-desktop/)

**Windows:**

Download and install [Docker Desktop for Windows](https://www.docker.com/products/docker-desktop/)

### Configure Docker Permissions (Linux)

Add your user to the docker group to run without sudo:

```
sudo usermod -aG docker $USER
newgrp docker
```

Log out and back in for changes to take effect.

### Verify Docker Installation

```
docker --version
docker run hello-world
```

## Available CHAMP Docker Images

CHAMP provides several pre-built images with different compiler and library configurations:

### GNU Compiler Images

| Image Tag | Description |
|-----------|-------------|
| `neelravi/champ:2.3.0` | CHAMP version 2.3.0 with GNU compilers |
| `neelravi/champ:gnu` | Latest CHAMP build with GNU compilers |
| `neelravi/champ:gnu-trexio` | GNU build with TREXIO support |

### Intel Compiler Images

| Image Tag | Description |
|-----------|-------------|
| `neelravi/champ:latest` | Latest CHAMP build with Intel oneAPI compilers |
| `neelravi/champ:intel` | CHAMP with Intel oneAPI compilers |
| `neelravi/champ:intel-trexio` | Intel build with TREXIO support |

## Pulling Docker Images

### Pull Specific Image

**Latest version with Intel compilers:**
```
docker pull neelravi/champ:latest
```

**GNU compiler version:**
```
docker pull neelravi/champ:gnu
```

**With TREXIO support (Intel):**
```
docker pull neelravi/champ:intel-trexio
```

**With TREXIO support (GNU):**
```
docker pull neelravi/champ:gnu-trexio
```

**Specific version (2.3.0):**
```
docker pull neelravi/champ:2.3.0
```

### List Downloaded Images

```
docker images | grep champ
```

## Running CHAMP in Docker

### Basic Usage

**Interactive shell:**
```
docker run -it neelravi/champ:latest /bin/bash
```

Once inside the container:
```
# Check CHAMP installation
which vmc.mov1
vmc.mov1 --version

# Navigate and run calculations
cd /path/to/workdir
vmc.mov1 -i input.inp -o output.out -e error
```

### Running Calculations with Volume Mounting

Mount your local directory to access input files and save output:

```
docker run -it -v /path/to/local/data:/data neelravi/champ:latest /bin/bash
```

Inside the container:
```
cd /data
vmc.mov1 -i input.inp -o output.out -e error
```

**Example:**
```
# Mount current directory
docker run -it -v $(pwd):/data neelravi/champ:gnu-trexio /bin/bash

# Inside container
cd /data
vmc.mov1 -i vmc.inp -o vmc.out -e error
```

### Running Single Command

Execute CHAMP directly without interactive shell:

```
docker run -v $(pwd):/data neelravi/champ:latest \
  vmc.mov1 -i /data/input.inp -o /data/output.out -e /data/error
```

### Using with MPI

For parallel calculations, use `mpirun` inside the container:

```
docker run -v $(pwd):/data neelravi/champ:intel \
  mpirun -np 8 vmc.mov1 -i /data/input.inp -o /data/output.out -e /data/error
```

## Common Workflows

### VMC Calculation

```
# Create a working directory
mkdir -p champ_calc
cd champ_calc

# Copy your input files here (vmc.inp, etc.)

# Run VMC
docker run -v $(pwd):/data neelravi/champ:gnu \
  vmc.mov1 -i /data/vmc.inp -o /data/vmc.out -e /data/error
```

### DMC Calculation

```
# Run VMC first to generate configurations
docker run -v $(pwd):/data neelravi/champ:latest \
  vmc.mov1 -i /data/vmc.inp -o /data/vmc.out -e /data/error

# Then run DMC
docker run -v $(pwd):/data neelravi/champ:latest \
  dmc.mov1 -i /data/dmc.inp -o /data/dmc.out -e /data/error
```

### Using TREXIO Files

For calculations using TREXIO format input:

```
docker run -v $(pwd):/data neelravi/champ:intel-trexio \
  vmc.mov1 -i /data/vmc_trexio.inp -o /data/vmc.out -e /data/error
```

## Advanced Docker Usage

### Creating Shell Alias

For easier usage, create an alias:

```
alias champ-vmc='docker run -v $(pwd):/data neelravi/champ:latest vmc.mov1'
alias champ-dmc='docker run -v $(pwd):/data neelravi/champ:latest dmc.mov1'
```

Add to `~/.bashrc` for persistence:
```
echo "alias champ-vmc='docker run -v \$(pwd):/data neelravi/champ:latest vmc.mov1'" >> ~/.bashrc
```

Then use simply:
```
champ-vmc -i /data/input.inp -o /data/output.out -e /data/error
```

### Inspecting Container Contents

```
# List files in container
docker run neelravi/champ:latest ls -lh /usr/local/bin

# Check library versions
docker run neelravi/champ:intel-trexio ldd /usr/local/bin/vmc.mov1

# Inspect environment
docker run neelravi/champ:latest env
```

### Saving Container State

If you make changes inside a container and want to save them:

```
# Start container
docker run -it --name my-champ neelravi/champ:latest /bin/bash

# Make changes inside...
# Exit container

# Commit changes
docker commit my-champ my-champ-custom:v1

# Run your custom image
docker run -it my-champ-custom:v1 /bin/bash
```

### Building Custom Image

You can create your own Docker image based on the CHAMP images or build from scratch.

#### Extending Existing CHAMP Image

Create a `Dockerfile` to add additional tools:

```dockerfile
FROM neelravi/champ:intel-trexio

# Install additional tools
RUN apt-get update && apt-get install -y \
    python3-pip \
    python3-numpy \
    python3-matplotlib

# Copy custom scripts
COPY my_scripts/ /usr/local/scripts/

# Set working directory
WORKDIR /data

CMD ["/bin/bash"]
```

Build and run:
```
docker build -t my-champ:custom .
docker run -it -v $(pwd):/data my-champ:custom
```

#### Building CHAMP from Scratch (GNU Compilers)

Here's the complete Dockerfile used to build the GNU-based CHAMP images:

```dockerfile
# Start from the latest Ubuntu image
FROM ubuntu:latest

# Set work directory
WORKDIR /app

# Install necessary packages
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    cmake \
    gnupg \
    autoconf \
    libtool \
    python3 \
    emacs \
    gcc \
    g++ \
    gfortran \
    openmpi-bin \
    libopenmpi-dev \
    gawk \
    build-essential \
    libscalapack-openmpi-dev \
    liblapack-dev \
    && rm -rf /var/lib/apt/lists/*

# Install HDF5 with MPI support
RUN curl -fsSLO 'https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads/manual/HDF5/HDF5_1_14_3/src/hdf5-1.14.3.tar' \
    && tar -xvf hdf5-1.14.3.tar \
    && cd hdf5-1.14.3 \
    && ./configure --prefix=/usr/local --enable-fortran --enable-parallel --enable-hl \
        FC=mpif90 CC=mpicc CXX=mpicxx \
    && make -j $(nproc) \
    && make install \
    && cd .. \
    && rm -rf hdf5-1.14.3 hdf5-1.14.3.tar

# Install TREXIO
RUN curl -fsSLO 'https://github.com/TREX-CoE/trexio/archive/refs/tags/v2.4.2.tar.gz' \
    && tar -xzvf v2.4.2.tar.gz \
    && cd trexio-2.4.2 \
    && cmake -S. -Bbuild \
        -DCMAKE_Fortran_COMPILER=mpif90 \
        -DCMAKE_C_COMPILER=mpicc \
    && cd build \
    && make -j $(nproc) \
    && make install \
    && cd ../.. \
    && rm -rf trexio-2.4.2 v2.4.2.tar.gz

# Install CHAMP
RUN git clone https://github.com/filippi-claudia/champ.git \
    && cd champ \
    && cmake -S. -Bbuild \
        -DCMAKE_Fortran_COMPILER=mpif90 \
        -DCMAKE_C_COMPILER=mpicc \
        -DENABLE_TREXIO=ON \
        -DTREXIO_INCLUDE_DIR=/usr/local/include \
        -DTREXIO_LIBRARY=/usr/local/lib/libtrexio.so \
    && cd build \
    && make -j $(nproc) \
    && cp /app/champ/bin/vmc.mov1 /app/champ/bin/dmc.mov1 /usr/local/bin/ \
    && cp /app/champ/bin/libparser.so /app/champ/bin/libpspline.a /usr/local/lib/

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
```

**Key components of the build:**

1. **Base image**: Ubuntu latest for package availability
2. **Build tools**: GCC, gfortran, OpenMPI, CMake
3. **Linear algebra**: OpenBLAS/LAPACK and ScaLAPACK
4. **HDF5**: Built with parallel support (--enable-parallel)
5. **TREXIO**: Version 2.4.2 with CMake build system
6. **CHAMP**: Built with TREXIO support enabled

**Build your own image:**

```
# Save Dockerfile
docker build -t my-champ:gnu-trexio .

# Run the image
docker run -it -v $(pwd):/data my-champ:gnu-trexio
```

**Notes:**
- The build process takes 20-30 minutes depending on your system
- The resulting image is approximately 2-3 GB
- All libraries are installed to `/usr/local`
- CHAMP executables are in `/usr/local/bin`

## Performance Considerations

### CPU Cores

By default, Docker uses all available CPU cores. To limit:

```
docker run --cpus=8 -v $(pwd):/data neelravi/champ:latest \
  mpirun -np 8 vmc.mov1 -i /data/input.inp -o /data/output.out
```

### Memory Limits

Set memory limit:

```
docker run -m 16g -v $(pwd):/data neelravi/champ:latest \
  vmc.mov1 -i /data/input.inp -o /data/output.out
```

### Temporary Files

Use tmpfs for temporary files:

```
docker run --tmpfs /tmp:rw,size=10g -v $(pwd):/data neelravi/champ:latest \
  vmc.mov1 -i /data/input.inp -o /data/output.out
```

## Updating Docker Images

Pull the latest version:

```
docker pull neelravi/champ:latest
```

Remove old images:

```
# List all CHAMP images
docker images | grep neelravi/champ

# Remove specific image
docker rmi neelravi/champ:old-tag

# Clean up unused images
docker image prune
```

## Docker on HPC Systems

Some HPC systems support Docker alternatives like Singularity/Apptainer that can run Docker images:

### Converting to Singularity

```
# Pull Docker image to Singularity format
singularity pull docker://neelravi/champ:latest

# Run with Singularity
singularity exec champ_latest.sif vmc.mov1 -i input.inp -o output.out
```

### Using in SLURM Jobs

```bash
#!/bin/bash
#SBATCH --job-name=champ-docker
#SBATCH --ntasks=16
#SBATCH --time=01:00:00

# Using Singularity on HPC
singularity exec champ_latest.sif \
  mpirun -np 16 vmc.mov1 -i input.inp -o output.out -e error
```

## Troubleshooting

### Permission Denied

If you get permission errors accessing files:

```
# Run with user ID mapping
docker run --user $(id -u):$(id -g) -v $(pwd):/data neelravi/champ:latest \
  vmc.mov1 -i /data/input.inp -o /data/output.out
```

### Container Exits Immediately

Check if command syntax is correct:

```
# Debug mode
docker run -it neelravi/champ:latest /bin/bash
# Then run commands manually
```

### Cannot Connect to Docker Daemon

```
# Start Docker service
sudo systemctl start docker

# Check status
sudo systemctl status docker
```

### Out of Disk Space

Clean up unused containers and images:

```
# Remove stopped containers
docker container prune

# Remove unused images
docker image prune -a

# Check disk usage
docker system df
```

## Advantages of Using Docker

- **Quick setup** - No compilation or dependency installation
- **Reproducibility** - Same environment across different machines
- **Version control** - Easy to switch between CHAMP versions
- **Portability** - Works on Linux, macOS, and Windows
- **Testing** - Try different configurations without affecting system
- **Clean removal** - Simply delete containers and images

## Limitations

- **Performance overhead** - Slight overhead compared to native execution
- **GPU support** - Requires additional configuration (NVIDIA Docker runtime)
- **HPC limitations** - Not all HPC systems support Docker (use Singularity instead)
- **Large image size** - Container images can be several GB

## Choosing the Right Image

| Use Case | Recommended Image |
|----------|-------------------|
| Quick testing | `neelravi/champ:latest` |
| GNU toolchain preference | `neelravi/champ:gnu` |
| Intel optimization | `neelravi/champ:intel` |
| TREXIO file input | `neelravi/champ:intel-trexio` or `gnu-trexio` |
| Specific version | `neelravi/champ:2.3.0` |
| Production calculations | `intel-trexio` (best performance + features) |

## Additional Resources

[:fontawesome-solid-book: Docker Documentation   :fontawesome-solid-arrow-up-right-from-square:](https://docs.docker.com/)

[:fontawesome-brands-docker: Docker Hub - CHAMP Images :fontawesome-solid-arrow-up-right-from-square:](https://hub.docker.com/r/neelravi/champ)

[:fontawesome-solid-book: Docker Best Practices :fontawesome-solid-arrow-up-right-from-square:](https://docs.docker.com/develop/dev-best-practices/)

[:fontawesome-solid-book: Singularity Documentation :fontawesome-solid-arrow-up-right-from-square:](https://sylabs.io/docs/)



## Next Steps

After pulling and testing a Docker image:

1. Verify the installation works with a simple test case
2. Learn about the [Command-Line Interface](supercomputers/cli.md)
3. Follow the [Tutorials](../../tutorials/index.md) for example calculations
4. Prepare your [Input Files](../../preparation/index.md)

For production calculations on HPC systems, consider [building from source](from_source.md) or using [Spack](using_spack.md)/[EasyBuild](using_easybuild.md) for optimal performance. 