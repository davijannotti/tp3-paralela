#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void generate_input(int N, int N_STEPS, double dt, const char* filename) {
    FILE* file = fopen(filename, "w");
    if (!file) {
        perror("Error opening file");
        exit(1);
    }

    fprintf(file, "%d\n", N);
    fprintf(file, "%d\n", N_STEPS);
    fprintf(file, "%lf\n", dt);

    srand(time(NULL));

    for (int i = 0; i < N; i++) {
        double mass = 1.0 + (rand() % 100) / 10.0; // Mass between 1.0 and 11.0
        double x = (rand() % 1000) - 500;          // x between -500 and 500
        double y = (rand() % 1000) - 500;          // y between -500 and 500
        double vx = (rand() % 10) - 5;             // vx between -5 and 5
        double vy = (rand() % 10) - 5;             // vy between -5 and 5
        
        fprintf(file, "%.2lf %.2lf %.2lf %.2lf %.2lf\n", mass, x, y, vx, vy);
    }

    fclose(file);
    printf("Generated %s with N=%d, Steps=%d, dt=%.4lf\n", filename, N, N_STEPS, dt);
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        printf("Usage: %s <N> <N_STEPS> <dt> <output_file>\n", argv[0]);
        return 1;
    }

    int N = atoi(argv[1]);
    int N_STEPS = atoi(argv[2]);
    double dt = atof(argv[3]);
    const char* filename = argv[4];

    generate_input(N, N_STEPS, dt, filename);

    return 0;
}
