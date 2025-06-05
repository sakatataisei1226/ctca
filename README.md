# MPIEMSES3D
3D Plasma Electromagnetic Particle MPI Parallel Simulator with the OhHelp Library

**MPIEMSES3D** is a high-performance simulation tool designed for 3D plasma electromagnetic particle-in-cell (PIC) simulations.
Using **MPI (Message Passing Interface)** for parallel computing, this simulator efficiently handles complex electromagnetic particle dynamics in large-scale simulations.
The integration with **the OhHelp library** further enhances this tool by adding dynamic load-balancing capabilities that optimize performance across distributed computing environments.

## Features

- **3D Plasma Electromagnetic Particle Simulations**: Simulate interactions between electromagnetic fields and particles within a 3D domain.

- **MPI Parallelism**: Supports distributed memory parallelism using MPI, allowing simulations to scale across multiple processors.

- **Dynamic Load Balancing**: Powered by the OhHelp library, `MPIEMSES3D` dynamically distributes computation across processors, improving performance in heterogeneous or unbalanced systems.

- **Support for Custom Boundary Conditions**: Incorporates boundary condition handling using another library `finbound`.
- **Coupled Computations**: Supports coupled computations with other programs using the CoToCoA framework.

## Clone
This repository has submodules. To clone the repository, use the `--recursive` option to clone the submodules as well.
Run the following command:

```bash
git clone --recursive https://github.com/Nkzono99/MPIEMSES3D.git
```

## Install
The simulator has three modes: **normal**, **cotocoa requester**, and **cotocoa worker**.
To build the simulator, simply run the following command in the root directory.
The executable file (`mpiemses3D`) will be generated in `MPIEMSES3D/bin` upon successful compilation:

In **normal** mode:
```bash
make
```
In **cotocoa requester** mode:
```bash
make requester
```
In **cotocoa worker** mode:
```bash
make worker
```

Sample code for **cotocoa requester** mode is available [here](https://github.com/kmr-gks/EMSES-CoToCoA-sample.git).

## Dependencies

The simulator relies on the following libraries:

### [OhHelp Library Package: Version 1.1.1 (2015/10/23)](http://www.para.media.kyoto-u.ac.jp/ohhelp/) <sup>1</sup>

The OhHelp library provides scalable, domain-decomposing dynamic load balancing essential for efficient particle-in-cell simulations across multiple processors.
It is particularly well-suited for cases where computational loads vary significantly across the simulation domain.

Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
                         (ACCMS, Kyoto University)

### [futils](https://github.com/Nkzono99/futils)
`futils` is a collection of utility functions.

Copyright (c) 2021 Nkzono99

### [finbound](https://github.com/Nkzono99/finbound)

`finbound` provides tools for managing boundary conditions within simulations, enhancing the flexibility of setting up customized simulation environments.

Copyright (c) 2021 Nkzono99

### [CoToCoA](https://github.com/tnanri/cotocoa)

`CoToCoA` is a framework to enable coupled computations such as multi-scale simulations, and forked by kmr-gks.

original:
https://github.com/tnanri/cotocoa

fork:
https://github.com/kmr-gks/cotocoa

Developers of CoToCoA:
- Takeshi Nanri (Kyushu Univ., Japan)
- Yuto Katoh (Tohoku Univ., Japan)
- Keiichiro Fukazawa (Kyoto Univ., Japan)
- Yohei Miyake (Kobe Univ., Japan)
- Kazuya Nakazawa (Kobe Univ., Japan)
- Jingde Zhow (Kyoto Univ., Japan)
- Youhei Sunada (Kobe Univ., Japan)
- Haichao Zhao (Kyoto Univ., Japan)

## References
[1] H. Nakashima, Y. Miyake, H. Usui and Y. Omura. OhHelp: A Scalable DomainDecomposing Dynamic Load Balancing for Particle-in-Cell Simulations. In Proc. Intl.
Conf. Supercomputing, pp. 90–99, June 2009. 4

This paper discusses the underlying architecture and performance benefits of the OhHelp library, which provides the load balancing features that mpiemses3d leverages for optimized particle-in-cell simulations.
