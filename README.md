# evolveX
Master thesis

# Overview

This repository contains an implementation of an algorithm used to predict **evolution environments** that can increase the production of a desired metabolite through adaptive laboratory evolution.

The workflow consists of three main steps:
1. Identify minimal reaction targets required to increase production of a metabolite.
2. Generate candidate environments from combinations of carbon and nitrogen sources.
3. Evaluate which environment produces the strongest growth coupling with the target reactions.

---

## Target Identification

The algorithm first identifies the smallest set of reactions whose flux must change in order to increase the flux of a desired metabolite.

A **trait matrix** is constructed and linear programming is used to determine an optimal production level (`v_opt`) where the flux of the desired metabolite is fixed.

A **mixed-integer linear programming (MILP)** problem is then solved to find the smallest number of reactions that must deviate from their wild-type flux ranges.

Binary variables `y_i` are used to indicate whether a reaction deviates from the wild-type range:

- `y_i = 0` → reaction flux remains within wild-type bounds  
- `y_i = 1` → reaction is allowed to deviate from wild-type bounds but must remain within global lower/upper bounds

The optimization minimizes the number of reactions where `y_i = 1`, meaning the minimal number of metabolic changes required to increase production.

---

## Generation of Growth Supporting Niches

Possible evolution environments are constructed from a **condition map (`condMap`)** containing selected carbon and nitrogen sources.

All combinations of these nutrients are generated to create candidate metabolic niches.

---

## Evolvex score

The goal is to determine how small the flux of each target reaction can be while the organism still grows optimally in that environment.

The values are summed to produce an overall score for each environment.

The environment with the **largest absolute score** is considered the most strongly growth-coupled environment.

---



