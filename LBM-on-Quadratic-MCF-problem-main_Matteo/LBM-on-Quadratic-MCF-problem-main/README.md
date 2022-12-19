# LBM on Quadratic MCF problem

Project for "Computational Mathematics for Learning and Data Analysis" course (a.y. 2021-2022) at the University of Pisa. The project consists in the implementation of a dual approach to solve a convex Quadratic Min-Cost Flow problem, where the Lagrangian Dual is solved by a level bundle method. More details can be found in the [report](https://github.com/Simoniuss/LBM-on-Quadratic-MCF-problem/blob/main/report/_CM__Project_report.pdf).

## Project tour

The project is structured as follows:
- `src`: contains the source code of the project.
    - `qcnd`: contains the source code to generate quadratic MCF problem.
    - `data`: contains the the .dmx and .qfc files of the generated problems and the .mat files of the relative problems used in Matlab.
- `report`: contains report and papers used for the mathematical formulation of the problem.

## Getting started
This is an example on how to run the project locally.

## Installation
1. Clone the repository.:
    ```
    git clone https://github.com/Simoniuss/LBM-on-Quadratic-MCF-problem.git
    ```
2. Compile the problem generator:
    ```
    cd src/qcnd
    make
    ```
3. Generate the problem:
    ```
    sh generateMCF.sh m rho k cf cg s
    ```
    More details about these parameters can be found in [readme](https://github.com/Simoniuss/LBM-on-Quadratic-MCF-problem/blob/main/src/qcnd/readme.txt). This script generate 3 files: .par, .dmx, .qfc, that are converted in a .mat using the python script convertertomat.py. For the converter you may have installed python 3 and the following packages: Numpy and Scipy.
    ```
    pip install numpy scipy
    ```
4. The data are stored in [data](https://github.com/Simoniuss/LBM-on-Quadratic-MCF-problem/blob/main/src/data).
5. Until now you can find a test problem 1000 edges in MCF_1000.mat