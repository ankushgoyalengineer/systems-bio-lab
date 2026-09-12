import matplotlib.pyplot as plt

def simulate_michael_menten(
    s_initial: float,
    e_initial: float,
    es_initial: 0.0,
    p_initial: 0.0,
    p_approx_initial: 0.0, 
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

    for i in range(50000):
        time_coordinates.append(current_time)
        s_coordinates.append(s_initial)
        e_coordinates.append(e_initial)
        es_coordinates.append(es_initial)
        p_coordinates .append(p_initial)
        p_approx_coordinates.append(p_approx_initial)
        ds = -(k1 * e_initial * s_initial * dt) + (k2 * es_initial * dt)
        de = -(k1 * e_initial * s_initial * dt) + (k2 * es_initial * dt) + (k3 * es_initial * dt)
        des = (k1 * e_initial * s_initial * dt) - (k2 * es_initial* dt) - (k3 * es_initial * dt)
        dp = (k3 * es_initial * dt)
        vmax = k3 * (e_initial + es_initial)
        v = (vmax * s_initial)/(km + s_initial)
        dp_approx = v * dt
        
        s_initial += ds
        e_initial += de
        es_initial += des
        p_initial += dp
        p_approx_initial += dp_approx
        current_time += dt

    return time_coordinates, s_coordinates, e_coordinates, es_coordinates, p_coordinates, p_approx_coordinates, s_initial, e_initial, es_initial, p_initial, p_approx_initial

def plot_result(time_coordinates, s_coordinates, e_coordinates, es_coordinates, p_coordinates, p_approx_coordinates, s_initial, e_initial, es_initial, p_initial, p_approx_initial):
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