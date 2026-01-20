# Contact Tracing with Incomplete Information (Julia)

This repository contains a multi-threaded **Agent-Based Model (ABM)** implemented in **Julia** to simulate infectious disease transmission under **imperfect contact tracing**.  
The model explicitly incorporates **missing information** in both infectors and contacts and evaluates how different omission mechanisms affect epidemic dynamics.

The code was developed for research purposes and is intended to reproduce and extend simulation results related to **Infector Omission (IO)** and **Contact Omission (CO)** scenarios.

---

## Model Overview

The simulation tracks individuals through detailed daily interactions across multiple social layers:

- Household
- School classroom / Workplace
- Friend gatherings
- Random encounters (residence community)

Disease progression, viral shedding, testing, quarantine, isolation, and contact tracing are explicitly modeled at the **individual level**.

### Contact Tracing Scenarios
The repository includes three main simulation scripts:

- **IO (Infector Omission)**  
  - Some confirmed infectors fail to be reported or traced.
- **SCO (Selective Contact Omission)**  
  - Only contacts from specific social layers (e.g., friends, community) may be missed.
- **UCO (Universal Contact Omission)**  
  - Contacts across all social layers may fail to be notified.

---

## Repository Structure

```text
.
├── ct_main_IO.jl      # Main simulation: Infector Omission (IO)
├── ct_main_SCO.jl     # Main simulation: Selective Contact Omission (SCO)
├── ct_main_UCO.jl     # Main simulation: Universal Contact Omission (UCO)
├── ct_ftn.jl          # Core functions (IBM engine, contact tracing, disease dynamics)
└── README.md
