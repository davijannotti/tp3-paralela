#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define G 6.67e-11

typedef struct {
    double mass;
    double x, y;
    double vx, vy;
} Particle;

typedef struct {
    double mass;
    double x, y;
} ParticleData;

void calculate_forces(Particle* local_particles, int local_N, ParticleData* all_particles, int total_N, double* fx, double* fy) {
    #pragma omp parallel for
    for (int i = 0; i < local_N; i++) {
        fx[i] = 0.0;
        fy[i] = 0.0;
        for (int j = 0; j < total_N; j++) {
            
            double dx = all_particles[j].x - local_particles[i].x;
            double dy = all_particles[j].y - local_particles[i].y;
            double dist_sq = dx*dx + dy*dy;
            
            if (dist_sq < 1e-10) continue;

            double dist = sqrt(dist_sq);
            double force = (G * local_particles[i].mass * all_particles[j].mass) / (dist * dist);
            
            fx[i] += force * (dx / dist);
            fy[i] += force * (dy / dist);
        }
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 3) {
        if (rank == 0) printf("Usage: %s <input_file> <output_file>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    int N, N_STEPS;
    double dt;
    Particle* all_initial_particles = NULL;

        if (rank == 0) {
        FILE* input = fopen(argv[1], "r");
        if (!input) {
            perror("Error opening input file");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        fscanf(input, "%d", &N);
        fscanf(input, "%d", &N_STEPS);
        fscanf(input, "%lf", &dt);

        all_initial_particles = (Particle*)malloc(N * sizeof(Particle));
        for (int i = 0; i < N; i++) {
            fscanf(input, "%lf %lf %lf %lf %lf", 
                   &all_initial_particles[i].mass, &all_initial_particles[i].x, &all_initial_particles[i].y, 
                   &all_initial_particles[i].vx, &all_initial_particles[i].vy);
        }
        fclose(input);
    }

    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N_STEPS, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    int local_N = N / size;
    int remainder = N % size;
    int* sendcounts = (int*)malloc(size * sizeof(int));
    int* displs = (int*)malloc(size * sizeof(int));
    
    int* sendcounts = (int*)malloc(size * sizeof(int));
    int* displs = (int*)malloc(size * sizeof(int));
    
    int sum = 0;
    for (int i = 0; i < size; i++) {
        sendcounts[i] = N / size;
        if (i < remainder) sendcounts[i]++;
        displs[i] = sum;
        sum += sendcounts[i];
    }
    local_N = sendcounts[rank];

    Particle* local_particles = (Particle*)malloc(local_N * sizeof(Particle));
    
    MPI_Datatype MPI_PARTICLE;
    MPI_Type_contiguous(5, MPI_DOUBLE, &MPI_PARTICLE);
    MPI_Type_commit(&MPI_PARTICLE);
    MPI_Scatterv(all_initial_particles, sendcounts, displs, MPI_PARTICLE, 
                 local_particles, local_N, MPI_PARTICLE, 
                 0, MPI_COMM_WORLD);

    if (rank == 0) free(all_initial_particles);

    double* fx = (double*)malloc(local_N * sizeof(double));
    double* fy = (double*)malloc(local_N * sizeof(double));

    MPI_Datatype MPI_PARTICLEDATA;
    MPI_Type_contiguous(3, MPI_DOUBLE, &MPI_PARTICLEDATA);
    
    MPI_Datatype MPI_PARTICLEDATA;
    MPI_Type_contiguous(3, MPI_DOUBLE, &MPI_PARTICLEDATA);
    MPI_Type_commit(&MPI_PARTICLEDATA);

    ParticleData* all_particle_data = (ParticleData*)malloc(N * sizeof(ParticleData));
    ParticleData* local_particle_data = (ParticleData*)malloc(local_N * sizeof(ParticleData));

    double start_time = MPI_Wtime();

    for (int step = 0; step < N_STEPS; step++) {
        for (int i = 0; i < local_N; i++) {
            local_particle_data[i].mass = local_particles[i].mass;
            local_particle_data[i].x = local_particles[i].x;
            local_particle_data[i].y = local_particles[i].y;
        }

        MPI_Allgatherv(local_particle_data, local_N, MPI_PARTICLEDATA,
                       all_particle_data, sendcounts, displs, MPI_PARTICLEDATA,
                       MPI_COMM_WORLD);

        calculate_forces(local_particles, local_N, all_particle_data, N, fx, fy);

        for (int i = 0; i < local_N; i++) {
            double ax = fx[i] / local_particles[i].mass;
            double ay = fy[i] / local_particles[i].mass;

            local_particles[i].vx += ax * dt;
            local_particles[i].vy += ay * dt;

            local_particles[i].x += local_particles[i].vx * dt;
            local_particles[i].y += local_particles[i].vy * dt;
        }
    }

    double end_time = MPI_Wtime();
    if (rank == 0) {
        printf("Simulation completed in %.4f seconds.\n", end_time - start_time);
    }

    Particle* final_particles = NULL;
    if (rank == 0) {
        final_particles = (Particle*)malloc(N * sizeof(Particle));
    }

    MPI_Gatherv(local_particles, local_N, MPI_PARTICLE,
                final_particles, sendcounts, displs, MPI_PARTICLE,
                0, MPI_COMM_WORLD);

    if (rank == 0) {
        FILE* output = fopen(argv[2], "w");
        if (!output) {
            perror("Error opening output file");
        } else {
            fprintf(output, "%d\n%d\n%lf\n", N, N_STEPS, dt);
            for (int i = 0; i < N; i++) {
                fprintf(output, "%.2lf %.2lf %.2lf %.2lf %.2lf\n", 
                        final_particles[i].mass, final_particles[i].x, final_particles[i].y, 
                        final_particles[i].vx, final_particles[i].vy);
            }
            fclose(output);
        }
        free(final_particles);
    }

    free(local_particles);
    free(fx);
    free(fy);
    free(all_particle_data);
    free(local_particle_data);
    free(sendcounts);
    free(displs);
    MPI_Type_free(&MPI_PARTICLE);
    MPI_Type_free(&MPI_PARTICLEDATA);

    MPI_Finalize();
    return 0;
}
