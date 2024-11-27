#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define L 5                // Edge length of the lattice
#define N (L * L)          // Total number of spins. creates [L x L] lattice
#define XNN 1              // Nearest neighbor in x-direction
#define YNN L              // Nearest neighbor in y-direction

int s[N];                  // Lattice of spins
double prob[5];            // Array for acceptance probabilities
double beta;               // Inverse temperature

// Random number generator for [0, 1)
double drandom() {
    return rand() / (RAND_MAX + 1.0);
}

// Initializing acceptance probabilities
// this is done and stored in an array once to make the program faster
//  so that the program doesn't have to calculate exponential everytime.
void initialize_probabilities() {
    for (int i = 2; i < 5; i += 2) {
        prob[i] = exp(-2 * beta * i);
    }
}

// Initialize lattice with random spins (+1 or -1)
void initialize_lattice() {
    for (int i = 0; i < N; i++) {
        s[i] = (rand() % 2 == 0) ? 1 : -1;
    }
}

// Perform a single Monte Carlo sweep over the lattice
void sweep() {
    for (int k = 0; k < N; k++) {
        int i = (int)(N * drandom());
        int nn, sum = 0, delta;

        // Calculate the sum of neighboring spins (helical boundary conditions)
        if ((nn = i + XNN) >= N) nn -= N; // Right neighbor
        sum += s[nn];
        if ((nn = i - XNN) < 0) nn += N;  // Left neighbor
        sum += s[nn];
        if ((nn = i + YNN) >= N) nn -= N; // Top neighbor
        sum += s[nn];
        if ((nn = i - YNN) < 0) nn += N;  // Bottom neighbor
        sum += s[nn];

        // Calculate the energy change
        delta = sum * s[i];

        // Decide whether to flip the spin
        // always flip the spin if Ev <= Eu
        // calculate the acceptance ratio for Ev > Eu
        if (delta <= 0 || drandom() < prob[delta]) {
            s[i] = -s[i];
        }
    }
}

// Calculate magnetization per spin
// magnetization = sum over all spins
// magnetization per spin = sum/N
double calculate_magnetization() {
    int sum = 0;
    for (int i = 0; i < N; i++) {
        sum += s[i];
    }
    return abs(sum) / (double)N;
}

// Calculate total energy of the system
double calculate_energy() {
    double energy = 0.0;
    for (int i = 0; i < N; i++) {
        int nn;
        // Right neighbor
        if ((nn = i + XNN) >= N) nn -= N;
        energy -= s[i] * s[nn];
        // Bottom neighbor
        if ((nn = i + YNN) >= N) nn -= N;
        energy -= s[i] * s[nn];
    }
    return energy;
}
// this only takes the right and left neighbors to ensure there is no repetition
// if we take all 4 neighbors, we simple return energy/2.0

// Main function
int main() {
    srand(time(NULL));  // Seed random number generator

    // Simulating between 0.1 - 5.0 temepratures
    double T_start = 0.1, T_end = 5.0, T_step = 0.1;
    int equilibration_steps = 2000;
    int measurement_steps = 18000;

    printf("T\t<m>\t\t<c>\n");
    for (double T = T_start; T <= T_end; T += T_step) {
        beta = 1.0 / T;
        initialize_lattice();
        initialize_probabilities();

        // Equilibrating the system
        for (int i = 0; i < equilibration_steps; i++) {
            sweep();
        }

        // Measuring properties
        double m_sum = 0.0, e_sum = 0.0, e2_sum = 0.0;
        for (int i = 0; i < measurement_steps; i++) {
            sweep();
            double e = calculate_energy();
            double m = calculate_magnetization();
            m_sum += m;
            e_sum += e;
            e2_sum += e * e;
        }

        // Calculate averages
        double m_avg = m_sum / measurement_steps;
        double e_avg = e_sum / measurement_steps;
        double e2_avg = e2_sum / measurement_steps;
        double c = (beta * beta / N) * (e2_avg - e_avg * e_avg);

        // Print results
        printf("%.2f\t%.4f\t\t%.4f\n", T, m_avg, c);
    }
    return 0;
}