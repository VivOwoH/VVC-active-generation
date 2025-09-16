# VVC-active-generation

Advanced voltage and reactive power control for distribution networks with high renewable energy penetration.

## Overview

Addresses voltage regulation challenges caused by distributed generation (DG) and renewable energy sources. Traditional OLTC and capacitor bank control is insufficient for managing rapid PV-induced voltage fluctuations.

## Structure

### Part 1: Day-Ahead Dispatch with Conventional VVC
- MISOCP optimization framework
- Three case studies:
  - Case 1: OLTC + CB (small load: 0.2-0.5 p.u.)
  - Case 2: OLTC + CB (medium load: 0.5-1.0 p.u.) 
  - Case 3: OLTC + CB + PV (medium load: 0.5-1.0 p.u.)

### Part 2: Enhanced Multi-Timescale VVC
- **Level 1**: Fast control (30s) - PV inverters, SVCs
- **Level 2**: Coordination (5min) - Voltage-based dispatch
- **Level 3**: Slow control (1h) - OLTC, capacitor banks

## System Model
- IEEE 33-bus distribution network
- OLTC: ±12.5% regulation, 10 tap positions
- Capacitor Banks: 5 × 900 kVA units
- PV Systems: 600 kVA total (~46% penetration)
- SVCs: 2 × ±900 kVA units (enhanced system)

## Key Results
- 60% reduction in CB reactive power requirements with PV integration
- Improved voltage regulation (0.95-1.05 p.u.) under high renewable penetration
- 30-second response times for voltage deviations

## Requirements
- MATLAB with Optimization Toolbox
- SOCP solver (MOSEK, Gurobi)
- IEEE 33-bus system data

## Applications
- Distribution network planning with renewables
- Real-time voltage control in smart grids
- Power quality improvement
- Grid modernization studies
