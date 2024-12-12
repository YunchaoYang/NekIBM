### 1. clone repo

```shell
git clone https://github.com/YunchaoYang/NekIBM.git
```
### 2.navigate to a sample case folder in IBM_test

```bash
cd NekIBM/IBM_test/simple
```

### 3. load relevant compiles

Note: A compatible list of compilers is not yet well tested. There meybe new updated version.

```bash
ml intel/2020.0.166   openmpi/4.1.6
```

### 4. compile the case
replace the uniform with the case name

```bash
./makenek uniform
```

### 5. modify uniform.batch and submit slurm jobs

Modify the uniform.batch file with the requested resources

```
sbatch uniform.batch
```

### 6.postprocess
The output particle phase files are saved in the vtk format, which could be opend by paraview or visit.
Run visnek to generate a meta file for the fluid phase. It can be loaded by paraview or visit.
