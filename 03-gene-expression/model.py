import matplotlib.pyplot as plt

class GeneNetwork:
    def __init__(self, n, alpha=10.0, kd=100.0, beta=5.0, ym=0.14, yp=0.02, dt=0.01, mrna_initial_conc=0.0, protein_initial_conc=0.0):
        self.alpha = alpha
        self.kd = kd
        self.n = n
        self.beta = beta
        self.ym = ym
        self.yp = yp
        self.dt = dt
        self.mrna_initial_conc = mrna_initial_conc
        self.protein_initial_conc = protein_initial_conc
        self.time_data = []
        self.mrna_data = []
        self.protein_data = []

    def simulate(self, steps):
        current_mrna_conc = self.mrna_initial_conc
        current_protein_conc = self.protein_initial_conc
        
        for i in range(steps):
            self.time_data.append(i * self.dt)
            self.mrna_data.append(current_mrna_conc)
            self.protein_data.append(current_protein_conc)
            
            if self.n == 0:
                transcription_rate = self.alpha
            else:
                transcription_rate = self.alpha / (1 + (current_protein_conc / self.kd) ** self.n)
            
            # Corrected variable names to reflect rates of change
            dm_dt = transcription_rate - (self.ym * current_mrna_conc)
            dp_dt = (self.beta * current_mrna_conc) - (self.yp * current_protein_conc)

            current_mrna_conc += dm_dt * self.dt
            current_protein_conc += dp_dt * self.dt

        return self.time_data, self.mrna_data, self.protein_data

    def plotter(self):
        plt.figure(figsize=(8, 5))
        plt.plot(self.time_data, self.mrna_data, color='blue', label='mRNA')
        plt.plot(self.time_data, self.protein_data, color='orange', label='Protein')
        plt.xlabel("Time")
        plt.ylabel("Concentration")
        plt.title(f"Gene Expression (Hill Coefficient n={self.n})")
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.savefig(f'03-gene-expression/plot_n{self.n}.png', dpi=150, bbox_inches='tight')
        plt.show()
        
def sweep_overshoot(n_values, steps=144000, **fixed_params):
    results = []
    for item in n_values:
        new_gene = GeneNetwork(n=item, **fixed_params)
        new_gene.simulate(steps)
        ratio = max(new_gene.protein_data) / (new_gene.protein_data)[-1]
        results.append((item, ratio))

    item_n, ratios = zip(*results)
    plt.figure(figsize=(8, 5))
    plt.plot(item_n, ratios, color='blue', marker='o')
    plt.xlabel('Hill coefficient (n)')
    plt.ylabel('Overshoot ratio')
    plt.title('Feedback Steepness vs Overshoot Ratio')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.show()
    return results

if __name__ == "__main__":
    # Define all standard parameters in one dictionary
    shared_parameters = {
        "alpha": 10.0,
        "kd": 100.0,
        "beta": 5.0,
        "ym": 0.14,
        "yp": 0.02,
        "dt": 0.01,
        "mrna_initial_conc": 0.0,
        "protein_initial_conc": 0.0
    }

    # 1. Run Baseline (No Feedback)
    print("Simulating baseline (n=0)")
    baseline_sim = GeneNetwork(n=0, **shared_parameters)
    baseline_sim.simulate(144000) 
    baseline_sim.plotter()

    # 2. Run the Damped Oscillator (n=4)
    print("Simulating damped oscillation (n=4)")
    damped_sim = GeneNetwork(n=4, **shared_parameters)
    damped_sim.simulate(144000)
    damped_sim.plotter()

    # 3. Run the Quantitative Sweep
    print("Running parameter sweep (n=1 to 10)")
    n_test_range = [1, 2, 4, 8, 10]
    overshooting = sweep_overshoot(n_values=n_test_range, steps=144000, **shared_parameters)
    print(overshooting)