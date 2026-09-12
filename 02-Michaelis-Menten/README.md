# Building a Mechanistic Enzyme Kinetics Engine

## Why I Built This
Coming from a life sciences background, I've spent a lot of time working around biological systems, but I learn best by actually building things from the ground up. Instead of just memorizing the standard Michaelis-Menten (M-M) equation for a portfolio project, I wanted to computationally derive it myself. 

The goal of this project was to build the raw, undisputed physics of a chemical reaction using Ordinary Differential Equations (ODEs) as a "ground truth", and then overlay the Michaelis-Menten approximation to physically see exactly when the math works, and more importantly, when it completely breaks down.

## The Process & The Math
I started by drafting the logic in a Jupyter Notebook so I could run the engine step-by-step before modularizing it into a final `.py` script. 

The reaction follows this standard path:
$$E + S \rightleftharpoons ES \rightarrow E + P$$

I built the engine around three specific rate constants: $k_1$ (binding), $k_2$ (dissociation), and $k_3$ (catalysis). 
One of the biggest early takeaways was strictly enforcing the Law of Mass Action and the conservation of mass. If an enzyme and substrate bind, they are *consumed* from the free pool (requiring negative signs in the ODEs), ensuring that the total enzyme in the system never magically increases ($E_{total} = E + ES$).

I broke the project into three distinct tasks:

### Task 1: The Ground Truth (Mechanistic Model)
I wrote a continuous-time simulation calculating the instantaneous rates of change for Substrate ($S$), Enzyme ($E$), Complex ($ES$), and Product ($P$). This required carefully tracking the negative and positive flow of molecules to prevent numerical leaks.

![Mechanistic Model Ground Truth](Mechanistic_model.png)


### Task 2: Proving the Approximation (Convergence)
To derive the Michaelis-Menten shortcut, I needed to prove the **Quasi-Steady-State Assumption (QSSA)**. By starting the simulation with $100$ Substrate and only $5$ Enzyme ($S \gg E$), the enzymes instantly saturated. 

As seen in the green line in the graph below, the $ES$ complex immediately spikes and holds a completely flat plateau ($\frac{d[ES]}{dt} = 0$). Because I mathematically proved this plateau exists, the black Michaelis-Menten approximation line perfectly tracks the grey ground-truth line. 

> ![Convergence: S >> E](MM_convergence.png)

### Task 3: Breaking the Math (Divergence)
To show exactly when standard biology textbook math fails, I intentionally sabotaged the initial conditions. I flooded the system so that Enzyme and Substrate concentrations were equal ($S = E$). 

Because the substrate didn't vastly outnumber the enzyme, the $ES$ complex never stabilized into a plateau. The QSSA condition was violated, and the M-M approximation (which is blind to the substrate hidden inside the $ES$ complex) incorrectly assumed the reaction was over, while the mechanistic ground-truth model continued to calculate the correct product formation.
> ![Divergence: S = E](MM_divergence.png)

## The Core Engine
Instead of relying on pre-built solvers, I built the numerical integration loop from scratch. Here is the core mathematical engine inside the simulation that calculates the instantaneous rates of change while strictly enforcing mass conservation. 

*(Full implementation and visualization logic are available in `model.py` and the Jupyter notebook).*

```python
        # 1. Mechanistic ODEs (Ground Truth)
        ds = -(k1 * e_initial * s_initial * dt) + (k2 * es_initial * dt)
        de = -(k1 * e_initial * s_initial * dt) + (k2 * es_initial * dt) + (k3 * es_initial * dt)
        des = (k1 * e_initial * s_initial * dt) - (k2 * es_initial * dt) - (k3 * es_initial * dt)
        dp = (k3 * es_initial * dt)
        
        # 2. Michaelis-Menten Approximation
        vmax = k3 * (e_initial + es_initial)
        v = (vmax * s_initial) / (km + s_initial)
        dp_approx = v * dt