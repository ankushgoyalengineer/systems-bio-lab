import matplotlib.pyplot as plt

def simulate_michael_menten(
    s_initial: float,
    e_initial: float,
    es_initial: float = 0.0,
    p_initial: float = 0.0,
    p_approx_initial: float = 0.0, 
    dt: float = 0.001, 
    k1: float = 10.0,
    k2: float = 2.0,
    k3: float = 2.0,
    steps: int = 10000
):
    print(f"s_initial={s_initial}, e_initial={e_initial}, k1={k1}, k2={k2}, k3={k3}")
    km = (k2 + k3) / k1
    print(f"km={km}")
    current_time = 0.0
    s_current = s_initial
    e_current = e_initial
    es_current = es_initial
    p_current = p_initial
    p_approx_current = p_approx_initial
    time_coordinates = []
    s_coordinates = []
    e_coordinates = []
    es_coordinates = []
    p_coordinates = []
    p_approx_coordinates = []
    print(f"e + es at start: {e_initial + es_initial}")

    for i in range(steps):
        time_coordinates.append(current_time)
        s_coordinates.append(s_current)
        e_coordinates.append(e_current)
        es_coordinates.append(es_current)
        p_coordinates .append(p_current)
        p_approx_coordinates.append(p_approx_current)
        ds = -(k1 * e_current * s_current * dt) + (k2 * es_current * dt)
        de = -(k1 * e_current * s_current * dt) + (k2 * es_current * dt) + (k3 * es_current * dt)
        des = (k1 * e_current * s_current * dt) - (k2 * es_current* dt) - (k3 * es_current * dt)
        dp = (k3 * es_current * dt)
        vmax = k3 * (e_current + es_current)
        v = (vmax * s_current)/(km + s_current)
        dp_approx = v * dt
        
        s_current += ds
        e_current += de
        es_current += des
        p_current += dp
        p_approx_current += dp_approx
        current_time += dt
    print(f"e + es at end: {e_current + es_current}")
    print(f"final product (mechanistic): {p_coordinates[-1]}")
    print(f"final product (approx): {p_approx_coordinates[-1]}")
    print(f"max substrate ever reached: {max(s_coordinates)}")
    return time_coordinates, s_coordinates, e_coordinates, es_coordinates, p_coordinates, p_approx_coordinates

def plot_result(time_coordinates, s_coordinates, e_coordinates, es_coordinates, p_coordinates, p_approx_coordinates):
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