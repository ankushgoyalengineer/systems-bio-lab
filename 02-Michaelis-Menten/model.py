import matplotlib.pyplot as plt

def simulate_michael_menten(
    s_initial: float,
    e_initial: float,
    es_initial: 0.0,
    p_initial: 0.0,
    p_approx: 0.0, 
    dt: float = 0.001, 
    k1: float = 10.0,
    k2: float = 2.0,
    k3: float = 2.0,
    steps: int = 10000
):
    km = (k2 + k3)/k1
    current_time = 0.0
    time_coordinates = []
    s_coordinates = []
    e_coordinates = []
    es_coordinates = []
    p_coordinates = []
    p_approx_coordinates = []

    for i in range(steps):
        time_coordinates.append(current_time)
        s_coordinates.append(s_conc)
        e_coordinates.append(e_conc)
        es_coordinates.append(es_conc)
        p_coordinates .append(p_conc)
        p_approx_coordinates.append(p_approx)
        ds = -(k1 * e_conc * s_conc * dt) + (k2 * es_conc * dt)
        de = -(k1 * e_conc * s_conc * dt) + (k2 * es_conc * dt) + (k3 * es_conc * dt)
        des = (k1 * e_conc * s_conc * dt) - (k2 * es_conc * dt) - (k3 * es_conc * dt)
        dp = (k3 * es_conc * dt)
        vmax = k3 * (e_conc + es_conc)
        v = (vmax * s_conc)/(km + s_conc)
        dp_approx = v * dt
        
        s_conc += ds
        e_conc += de
        es_conc += des
        p_conc += dp
        p_approx += dp_approx
        current_time += dt

    return time_coordinates, s_coordinates, e_coordinates, es_coordinates, p_coordinates, p_approx_coordinates, s_initial, e_initial, es_initial, p_initial

def plot_result(time_coordinates, s_coordinates, e_coordinates, es_coordinates, p_coordinates, p_approx_coordinates, s_initial, e_initial, es_initial, p_initial):
    plt.plot(time_coordinates, s_coordinates, label="Substrate concentration", color="blue")
    plt.plot(time_coordinates, e_coordinates, label="Enzyme concentration", color="red")
    plt.plot(time_coordinates, es_coordinates, label="ES Complex concentration", color="green")
    plt.plot(time_coordinates, p_coordinates, label="Product concentration", color="grey")
    plt.plot(time_coordinates, p_approx_coordinates, label="Product approx concentration", color="black")

    plt.xlabel("Time")
    plt.ylabel("Concentration (µM)")
    plt.title("Enzyme Kinetics Over Time")
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.show()