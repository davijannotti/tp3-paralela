# N-Body Simulation (MPI + OpenMP)

Este pacote contém a implementação de uma simulação de N-Corpos utilizando uma abordagem híbrida de MPI e OpenMP.

## Arquivos

- `nbody_seq.c`: Versão sequencial da simulação.
- `nbody_par.c`: Versão paralela (MPI + OpenMP).
- `input_gen.c`: Gerador de arquivos de entrada.
- `relatorio.pdf`: Relatório técnico do trabalho.

## Compilação

Para compilar os códigos, utilize os seguintes comandos:

```bash
# Gerador de entrada
gcc -o input_gen input_gen.c -lm

# Versão Sequencial
gcc -o nbody_seq nbody_seq.c -lm

# Versão Paralela (Requer biblioteca MPI instalada)
mpicc -o nbody_par nbody_par.c -lm -fopenmp
```

## Execução

### 1. Gerar Entrada
Primeiro, gere um arquivo de entrada com o número desejado de corpos (N), passos de tempo (STEPS) e dt:

```bash
# Uso: ./input_gen <N> <STEPS> <dt> <output_file>
./input_gen 1000 100 0.1 input.txt
```

### 2. Executar Sequencial

```bash
# Uso: ./nbody_seq <input_file> <output_file>
./nbody_seq input.txt output_seq.txt
```

### 3. Executar Paralelo

Para rodar a versão paralela, utilize o `mpirun` especificando o número de processos MPI (`-np`). O número de threads OpenMP por processo pode ser controlado pela variável de ambiente `OMP_NUM_THREADS`.

```bash
# Exemplo: 4 Processos MPI, cada um com 2 Threads OpenMP (Total 8 cores)
export OMP_NUM_THREADS=2
mpirun -np 4 ./nbody_par input.txt output_par.txt
```

## Verificação

Os arquivos de saída contêm o estado final das partículas. Você pode comparar `output_seq.txt` e `output_par.txt` para verificar a correção (pequenas diferenças de ponto flutuante são esperadas).
