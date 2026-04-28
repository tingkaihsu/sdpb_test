# Setup of SDPB problem in a new computer

## NTU CPU

1. log in command to NTU CPU
```bash

```
2. run `pmp2adp`
```bash
    docker run --shm-size=16384m -v "$(pwd):/data" bootstrapcollaboration/sdpb:master pmp2sdp 2048 -i /data/n_pmp.json -o /data/out
```

3. run `sdpb`
```bash
    docker run --shm-size=16384m -v "$(pwd):/data" bootstrapcollaboration/sdpb:master mpirun --allow-run-as-root -n 24 sdpb --writeSolution="x,y,z,X,Y" --maxComplementarity=1e1000 --dualityGapThreshold=1e-30 --stepLengthReduction=0.7 --primalErrorThreshold=1e-30 --dualErrorThreshold=1e-30 --precision=2048 --procsPerNode=32 --maxIterations=50000 -s /data/out
```