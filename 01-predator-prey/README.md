# Ecosystem Dynamics: Lotka-Volterra Predator-Prey Simulation

## Overview
This project is a computational biology simulation that models the continuous dynamic interactions between predator and prey populations. Built with Python, the simulation utilizes numerical integration (Euler's method) to solve coupled differential equations based on the classic Lotka-Volterra model. 

To reflect realistic physical environments, the base model is augmented with a logistic growth term representing the environmental carrying capacity ($K$). This introduces "environmental friction," allowing the simulation to demonstrate the mathematical transition from volatile, endless population cycles to damped oscillations and stable equilibria.

## Key Engineering Features
* **Simultaneous State Updating:** The engine computes instantaneous rates of change ($dx$ and $dy$) simultaneously before applying time-step updates, successfully preventing the numerical instability and overshoot errors common in continuous-time simulations.
* **Modular Architecture:** The codebase strictly separates the mathematical simulation engine from the data visualization layer, allowing for highly reusable code.
* **Dynamic Parameterization:** Ecosystem variables (growth rates, death rates, hunting efficiency, and carrying capacity) are exposed as configurable parameters, enabling rapid scenario testing and comparative analysis.

## Visual Artifacts & Analysis

### 1. The Impact of Carrying Capacity (Resource Limits)
By plotting two different environmental scenarios on the same time axis, we can observe how carrying capacity ($K$) acts as a physical failsafe. 
* **Orange (Unlimited Resources, K=1000):** Results in aggressive, continuous population swings.
* **Blue (Limited Resources, K=25):** The friction of a lower carrying capacity forces a damped oscillation, smoothing the population curves into a permanent, flat equilibrium.

> ![alt text](image-2.png)

### 2. Phase Portraits: Orbital Cycles vs. Equilibrium Spirals
By removing the time axis and plotting the prey population directly against the predator population, we can visualize the ecosystem's trajectory. 

With unlimited resources, the populations remain locked in a massive, continuous orbit. However, introducing a strict resource limit ($K=25$) causes the phase portrait to collapse into an inward spiral, pulling the populations directly into a mathematical equilibrium of exactly 10 prey and 6 predators.

> ![alt text](image-1.png)
> ![alt text](image.png)

## Tech Stack
* **Language:** Python 3.x
* **Libraries:** Matplotlib (Data Visualization)
* **Concepts:** Numerical Integration, Coupled Differential Equations, Data Visualization, Object-Oriented/Functional Design

## How to Run Locally

1. Clone this repository to your local machine:
   ```bash
   git clone [https://github.com/YourUsername/YourRepositoryName.git](https://github.com/YourUsername/YourRepositoryName.git)