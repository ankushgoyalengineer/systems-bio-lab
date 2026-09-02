import matplotlib.pyplot as plt

def simulate_lotka_volterra(
    prey_initial: float, 
    predator_initial: float, 
    dt: float = 0.001, 
    c: float = 1.0, 
    delta: float = 0.1, 
    alpha: float = 1.0, 
    beta: float = 0.1, 
    K: float = 25.0,
    steps: int = 10000
):
    """
    Simulates predator-prey dynamics using the Lotka-Volterra equations with logistic growth.
    
    Args:
        prey_initial: Starting prey population.
        predator_initial: Starting predator population.
        dt: Time step for numerical integration.
        c: Predator natural death rate.
        delta: Predator growth efficiency from eating prey.
        alpha: Prey natural growth rate.
        beta: Predation rate (prey hunted per predator).
        K: Environmental carrying capacity for prey.
        steps: Number of iterations for the simulation.
    """
    
    # Initialize tracking lists
    time_data = []
    prey_data = []
    predator_data = []
    
    current_prey = prey_initial
    current_predator = predator_initial
    current_time = 0.0

    # Single simulation loop handling all scenarios natively via the math
    for i in range(steps):
        time_data.append(current_time)
        prey_data.append(current_prey)
        predator_data.append(current_predator)

        # Calculate rates of change
        dy = (delta * current_prey * current_predator * dt) - (c * current_predator * dt)
        dx = (alpha * current_prey * (1 - (current_prey / K)) * dt) - (beta * current_prey * current_predator * dt)

        # Update populations
        current_predator += dy
        current_prey += dx
        current_time += dt

        # Prevent physically impossible negative populations
        if current_predator <= 0:
            current_predator = 0
        if current_prey <= 0:
            current_prey = 0


    return time_data, prey_data, predator_data, predator_initial, prey_initial

def plot_result(time_data, prey_data, predator_data, predator_initial, prey_initial):
    # Determine dynamic title based on starting conditions
    if predator_initial == 0:
        graph_title = "Growth of Prey Population (No Predators)"
    elif prey_initial == 0:
        graph_title = "Growth of Predator Population (No Prey)"
    else:
        graph_title = "Predator vs Prey Populations Over Time"    
    # Plotting the results
    plt.plot(time_data, prey_data, label="Prey Population", color="blue")
    plt.plot(time_data, predator_data, label="Predator Population", color="red")

    plt.xlabel("Time")
    plt.ylabel("Population Size")
    plt.title(graph_title)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6) # Added a light grid for visual polish
    plt.show()
