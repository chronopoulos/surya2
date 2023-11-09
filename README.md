# surya2
a parallelized version of surya, the mean-field solar dynamo simulation

## Dependencies
```
sudo apt-get install gfortran openmpi-bin libopenmpi3 libopenmpi-dev
```

## Build
```
cd parallel
make
```

## Run
```
mpirun -n $(nproc) surya
```
