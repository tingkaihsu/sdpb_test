# Usage
```bash
    sudo docker run -v /my/project/:/usr/local/share/sdpb bootstrapcollaboration/sdpb:master mpirun --allow-run-as-root -n 4 pmp2sdp --precision 1024 -i /usr/local/share/sdpb/pmp.json -o /usr/local/share/sdpb/sdp
```
This would convert a *polynomial matrix program* problem into an internal SDP format that SDPB can immediately solve.
```bash 
    docker run -v /my/project/:/usr/local/share/sdpb bootstrapcollaboration/sdpb:master mpirun --allow-run-as-root -n 4 sdpb --precision 1024 -s /usr/local/share/sdpb/sdp -o /usr/local/share/sdpb/sdpb
```
One can find the output of this problem in `sdpb/sdpb/out.txt`, and it shows the reason for termination, the final primal/dual objective value, the final duality gap, the final primal/dual errors, and the total runtime.
The output one may be interested in is the vector **y**, which is stored in `sdpb/sdpb/y.txt`.
