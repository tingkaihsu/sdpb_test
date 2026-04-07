# Discretization (Adaptive Refinement)
SDP problem:
$$
    \text{Subjected to the bound: } f(x) \ge 0
$$
$$
    \text{Maximize the objection: } \vec{b}(x) \cdot \vec{y} \text{where the normalization is } \vec{n}\cdot\vec{y} = 1
$$
The function $f(x)$ is given and one has to design the objection and the normalization according to the problem at hand.
We initially take some sampling points $x = [\epsilon_x, \epsilon_x+\delta_x, ..., 1-\delta_x]$. One straightforward consequence of only imposing positivity bound at a discrete set of $x$ is that the function with the solution is negative between two discretized values where the inequalities are saturated. 
*If the left-hand side of the inequality is zero at some values of $x$, it will generically be negative on one side of that zero.*

One can mitigate this problem by adaptively refining the discretization. Notice that the resulting functional would be negative between pairs of points in the initial discretization. **We identify these negative regions and add new positivity constraints in the finer-spaced grid that covers the negative regions.**

My collaborator suggests that *要跑完一次sdp 之後精進sampling points*, *你要回來確認你的functional是不是大於零*. His approach is that *另外寫一個shell script執行全部的東西, 包括讀檔分析*.

## My Suggestion
To understand my collaborator's approach, one has to study how SDPB solver works.
### PMP to SDP
First, we do NOT write the SDP problem directly. We start with writing down the PMP problem, and use SDPB solver to translate it into a SDP problem. By sampling, we mean the sampling points of the **PMP** problem, not the SDP problem. As we discussed before, for non-polynomial problem, we treat its numeric values at every sampling points as constant functions. For $3$ non-polynomial constraints with $10$ sampling points, we translate it into $3*10 = 30$ constant-function constraints at $10$ sampling points. 
### Next Step
The sampling points is taken in a `mathematica` code.
With the help of `SDPB.m`, the PMP problem written in the `mathematica` code output a `JSON` file with a certain format.

Then we copy the `JSON` file (`input.json`) to the path `/my/project/` with `sudo` permission, where we run the SDPB solver there. Here is the command
```bash
    sudo docker run -v /my/project/:/usr/local/share/sdpb bootstrapcollaboration/sdpb:master mpirun --allow-run-as-root -n 4 pmp2sdp --precision 1024 -i /usr/local/share/sdpb/input.json -o /usr/local/share/sdpb/sdp
```
This will generate a folder (`/my/project/sdp/`). 
IMPROTANT: One should make sure that there is no `/my/project/sdp/` of the previous run, DELETE IT otherwise the `bash` command would generate errors.

### SDPB solver
Finally, we successfully translate the PMP problem into the SDP problem to be solved. We can run the SDPB solver at the folder `/my/project` with `sudo` permission:
```bash
    sudo docker run -v /my/project/:/usr/local/share/sdpb bootstrapcollaboration/sdpb:master mpirun --allow-run-as-root -n 4 sdpb --precision=1024 -s /usr/local/share/sdpb/sdp
```
This will generate two folders: one for the checkpoints during the finding of optimizations, and the other for the record of solution (might fail due to some reasons, such as memory limit or time limit, etc.). Now 