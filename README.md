# Communication Topology Co-Design for Islanded Microgrids

MATLAB-based stability analysis, dynamic validation, and communication link-and-weight optimization for enhanced droop-free control of an islanded four-DER microgrid.

**Student:** Eileen Figueroa  
**Advisor:** Dr. Lei Wu  
**Program:** M.S. Electrical Engineering, Stevens Institute of Technology  
**Course:** EE 800 Engineering Project, Fall 2026

## Project Overview

Distributed energy resources (DERs) in an islanded microgrid must coordinate their active-power contributions while maintaining stable operation. Consensus-based droop-free controllers use a communication network to exchange information among DERs. Although a dense communication network can improve stability, it also increases communication cost and complexity.

This project evaluates original normalized active-power control (O-NAPC) and enhanced normalized active-power control (A-NAPC), studies how communication topology affects stability and dynamic performance, and identifies a sparse weighted communication network that satisfies specified engineering requirements.

## Research Question

How can the communication links and their weights be optimized to reduce communication requirements while preserving the stability and dynamic performance of an A-NAPC-controlled islanded microgrid?

## System Model

The study uses a four-DER islanded microgrid with fixed electrical-network parameters, controller gain, and DER capacities. Communication topology is represented by the weighted Laplacian matrix `L(w)`, where a zero edge weight indicates that the corresponding communication link is absent.

For the A-NAPC model, the closed-loop state matrix is

<p align="center"><strong>A(<i>w</i>) = −hB<sub>el</sub>D<sub>p</sub>L(<i>w</i>)D<sub>p</sub></strong></p>

The effective stability margin is determined from the nonzero closed-loop eigenvalues:

<p align="center"><strong>ρ(<i>w</i>) = −max Re{λ<sub>k</sub>(A(<i>w</i>))}</strong></p>

where the maximum is evaluated over the nonzero eigenvalues satisfying λ<sub>k</sub>(A(<i>w</i>)) ≠ 0.

## Methodology

The project progressed through twelve MATLAB experiments:

| Experiment | Purpose |
|---|---|
| 1 | Establish the four-DER A-NAPC baseline and calculate eigenvalues, stability margin, consensus error, and settling time |
| 2–4 | Compare O-NAPC and A-NAPC across communication topologies and controller gains |
| 5 | Evaluate sensitivity to asymmetric communication weights |
| 6 | Enumerate all 38 connected labeled graphs for the four-node system |
| 7 | Compare dynamic responses under common and localized DER 4 disturbances |
| 8 | Perform an initial communication cost–performance co-design study |
| 9 | Analyze modal excitation to explain disturbance-dependent settling behavior |
| 10 | Compare weighted-sum and constrained topology-selection formulations |
| 11 | Validate per-node settling thresholds using a maximum-node criterion |
| 12 | Jointly optimize communication-link selection and continuous edge weights |

The final optimization minimizes a smooth approximation of the number of active links and the total communication weight while enforcing:

- At least 98% of the complete-graph stability margin
- Common-input settling time no greater than 0.40 s
- Localized DER 4 settling time no greater than 0.45 s
- A connected communication graph

## Final Optimized Design

The final design contains four active communication links:

| Communication Link | Optimized Weight |
|---|---:|
| DER 1–DER 2 | 1.000 |
| DER 1–DER 3 | 1.000 |
| DER 1–DER 4 | 1.000 |
| DER 2–DER 4 | 0.562 |
| DER 2–DER 3 | 0 |
| DER 3–DER 4 | 0 |

### Comparison with the Unit-Weight Design

| Metric | Experiment 10 Unit Weights | Experiment 12 Optimized Weights |
|---|---:|---:|
| Active links | 4 | 4 |
| Total communication weight | 4.000 | 3.562 |
| Stability margin | 11.121 | 11.054 |
| Complete-graph margin retained | 98.597% | 98.000% |
| Common-input settling time | 0.238 s | 0.238 s |
| Localized DER 4 settling time | 0.187 s | 0.185 s |
| Algebraic connectivity | 1.000 | 1.000 |

The optimized design reduces total communication weight by approximately **10.9%** while satisfying all stability, dynamic-performance, and connectivity constraints.

## Repository Structure

```text
microgrid-communication-codesign/
├── README.md
├── .gitignore
└── src/
    ├── README.md
    ├── baseline/
    │   ├── README.md
    │   ├── asym_lattice.m
    │   ├── simu_4.m
    │   └── stabilitymargin.m
    └── experiments/
        ├── README.md
        ├── experiment1_simu_4_analysis.m
        ├── simu_4_controlled_experiments.m
        ├── experiment5_asymmetry_sensitivity.m
        ├── experiment6_topology_enumeration.m
        ├── experiment7_dynamic_scenarios.m
        ├── experiment8_codesign.m
        ├── Experiment9_Modal_Excitation_Analysis.m
        ├── Experiment10_Formal_Topology_Codesign.m
        ├── Experiment11_Per_Node_Threshold_Validation.m
        └── Experiment12_Joint_Link_Weight_Optimization.m
```

## MATLAB Requirements

The scripts use the following MATLAB functionality:

- MATLAB
- Control System Toolbox
- Optimization Toolbox for `fmincon`
- Statistics and Machine Learning Toolbox for `unifrnd` in the original asymmetric-lattice script
- Global Optimization Toolbox is optional; Experiment 12 can use `patternsearch` when available and otherwise falls back to multistart `fmincon`

## Running the Project

1. Clone or download the repository.
2. Open MATLAB.
3. Add `src/baseline` and `src/experiments` to the MATLAB path.
4. Run the experiment scripts individually in numerical order.
5. Run `Experiment12_Joint_Link_Weight_Optimization.m` to reproduce the final optimized communication design.

Experiment scripts save their figures as PNG files and numerical results as CSV files in the current MATLAB working directory.

## Attribution

The original `asym_lattice.m`, `simu_4.m`, and `stabilitymargin.m` baseline scripts were provided by Dr. Lei Wu. The experiment scripts extend those starting models through additional stability calculations, topology comparisons, disturbance analysis, design-space evaluation, constraint validation, and communication link-and-weight optimization.

## Primary References

1. W. Liu et al., “Stability Analysis of Droop-Free Controlled Islanded Microgrids with Asymmetric Communication Networks,” *Cyber-Physical Energy Systems*, 2026.

2. W. Liu et al., “On Theoretical Stability Proof and Stability Margin Analysis of Enhanced Droop-Free Control Schemes for Islanded Microgrids,” *IEEE Transactions on Industrial Informatics*, vol. 22, no. 6, 2026.

3. M. J. Najafirad and S. Welikala, “Dissipativity-Based Distributed Control and Communication Topology Co-Design for Nonlinear DC Microgrids,” arXiv:2603.27427, 2026.

## Project Status

The optimization and MATLAB validation stages are complete. Final reporting and research-poster preparation are in progress.
