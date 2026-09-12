# Particle Diffusion & Energy Dynamics Simulation

**Overview:**
A computational physics model simulating particle diffusion and energy configuration across an L x L lattice grid over time (represented as steps). This simulation tracks the decay of particle populations and calculates the average energy per particle, modeling how particles interact, annihilate, and reach equilibrium states based on bond energy limits and temperature values.

**Technical Stack:**
*   **Language:** Python
*   **Scientific Libraries:** NumPy, Matplotlib, Math, Random[cite: 2]
*   **Deployment:** UVA Rivanna HPC Cluster (SLURM workload manager)

**Core Simulation Mechanics**
*   **HPC Integration:** Utilizes `os` to fetch the `SLURM_JOB_ID` for cluster deployment and applies `matplotlib.use('Agg')` for headless data visualization without a graphical interface.
*   **Lattice Generation:** Initializes an L x L zero-matrix, randomly distributing an initial particle population based on a set probability density.
*   **Kinematic Logic:** Particles evaluate up to 6 directional moves per time step. Moves are evaluated using a boundary-aware neighbor-counting algorithm to determine bond energies between adjacent sites.
*   **Thermodynamic Acceptance:** Transition probabilities are calculated by comparing the energy difference between current and proposed states, utilizing a threshold criteria of A = min(1, e^(-E/T)).

**Generated Outputs:**
The script runs hundreds of parallel trials to generate time-series visualizations, including:
*   Individual trial staircase patterns for total configuration energy and particle counts.
*   Logarithmic plots of average particle decay compared against the theoretical 1/t expectation.
*   Average energy per particle tracking over defined time steps.
