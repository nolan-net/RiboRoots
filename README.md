# RiboRoots

## Overview

RiboRoots is a streamlined bioinformatics application focused on analyzing and visualizing the phylogenetic relationships of bacterial strains. This project integrates Python-based algorithms for the comparison of 16S rRNA gene sequences to construct detailed phylogenetic trees. Utilizing 'rrnDB version 5.8', a comprehensive dataset from the Center for Microbial Systems at the University of Michigan, RiboRoots includes data on several thousand bacterial strains and their corresponding 16S rRNA genes.

## Objective

The aim of RiboRoots is to develop a program using Python that allows users to select bacterial strains, compare their RNA sequences, and visually represent the relationships through phylogenetic trees.

## Key Features

1. **Reliable Data Source**: Incorporates 'rrnDB version 5.8', a dataset from the Center for Microbial Systems, University of Michigan.
2. **Algorithmic Sequence Analysis**: Uses Python to compare 16S rRNA gene sequences accurately.
3. **Command Line Interface**: Utilizes an easy to use CLI for use across various devices/hardware.
4. **Inclusive Design**: Ensures the application is intuitive and accessible for users with various levels of expertise.
5. **Docker Build**: Has a readily available Docker Image on Docker Hub (nolannet/riboroots)

## Setup
Currently the only way to use local graphing with ete3 package is to clones this repo and use the source code. However, an easy to use docker image is availabe to users less familiar with such a process.

### Docker
 * Image availabe here:
  - https://hub.docker.com/r/nolannet/riboroots

### Local Installation
 * You can either clone this repo, or just download the main.py file. This program is configured to install all necessary packages at user discretion which will have it all ready to go as long as a python environment is setup.

### Examples:
#### Intro to program + alignment
<img width="788" alt="Screenshot 2024-05-11 at 3 42 20 AM" src="https://github.com/nolan-net/RiboRoots/assets/91292872/70ffa27e-d927-4fc2-87a2-480b3f2d4e00">

#### Local Graph Using ete3
<img width="545" alt="Screenshot 2024-05-11 at 3 43 33 AM" src="https://github.com/nolan-net/RiboRoots/assets/91292872/7b6b6fba-3b7d-4e76-acfc-005b459edf17">

#### Graph of Legionella Strands Using Interactive Tree Of Life (separate graphing website: https://itol.embl.de/ )
<img width="755" alt="Screenshot 2024-05-11 at 3 44 58 AM" src="https://github.com/nolan-net/RiboRoots/assets/91292872/65f2a870-c720-48e8-9ee8-da791cb85a96">

## License

This project is open for general use under the GNU GENERAL PUBLIC LICENSE.

---

RiboRoots offers a practical and efficient solution for exploring and understanding the phylogenetic relationships among bacterial species, making it easier for researchers and enthusiasts to delve into the world of microbial phylogenetics.
