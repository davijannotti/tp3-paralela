#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define G 6.67e-11

typedef struct {
    double mass;
    double x, y;
    double vx, vy;
    double fx, fy;
} Particle;

void calculate_forces(Particle* particles, int N) {
    for (int i = 0; i < N; i++) {
        particles[i].fx = 0.0;
        particles[i].fy = 0.0;
        for (int j = 0; j < N; j++) {
            if (i == j) continue;

            double dx = particles[j].x - particles[i].x;
            double dy = particles[j].y - particles[i].y;
            double dist_sq = dx*dx + dy*dy;
            double dist = sqrt(dist_sq);
            
            if (dist < 1e-5) dist = 1e-5;

            double force = (G * particles[i].mass * particles[j].mass) / (dist * dist);
            
            particles[i].fx += force * (dx / dist);
            particles[i].fy += force * (dy / dist);
        }
    }
}

void update_particles(Particle* particles, int N, double dt) {
    for (int i = 0; i < N; i++) {
        double ax = particles[i].fx / particles[i].mass;
        double ay = particles[i].fy / particles[i].mass;

        particles[i].vx += ax * dt;
        particles[i].vy += ay * dt;

        particles[i].x += particles[i].vx * dt;
        particles[i].y += particles[i].vy * dt;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        printf("Usage: %s <input_file> <output_file>\n", argv[0]);
        return 1;
    }

    FILE* input = fopen(argv[1], "r");
    if (!input) {
        perror("Error opening input file");
        return 1;
    }

    int N, N_STEPS;
    double dt;
    fscanf(input, "%d", &N);
    fscanf(input, "%d", &N_STEPS);
    fscanf(input, "%lf", &dt);

    Particle* particles = (Particle*)malloc(N * sizeof(Particle));
    for (int i = 0; i < N; i++) {
        fscanf(input, "%lf %lf %lf %lf %lf", 
               &particles[i].mass, &particles[i].x, &particles[i].y, 
               &particles[i].vx, &particles[i].vy);
    }
    fclose(input);

    clock_t start = clock();

    for (int step = 0; step < N_STEPS; step++) {
        calculate_forces(particles, N);
        update_particles(particles, N, dt);
    }

    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Simulation completed in %.4f seconds.\n", time_spent);

    FILE* output = fopen(argv[2], "w");
    if (!output) {
        perror("Error opening output file");
        free(particles);
        return 1;
    }

    fprintf(output, "%d\n%d\n%lf\n", N, N_STEPS, dt);
    for (int i = 0; i < N; i++) {
        fprintf(output, "%.2lf %.2lf %.2lf %.2lf %.2lf\n", 
                particles[i].mass, particles[i].x, particles[i].y, 
                particles[i].vx, particles[i].vy);
    }
    fclose(output);
    free(particles);

    return 0;
}
