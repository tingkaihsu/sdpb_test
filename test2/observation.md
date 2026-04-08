# Discretization (Adaptive Refinement)
SDP problem:
$$
    \text{Subjected to the bound: } \vec{f}(x) \cdot \vec{y} \ge 0,
$$

$$
    \text{Maximize the objection: } \vec{b}(x) \cdot \vec{y} \text{, where the normalization is } \vec{n}\cdot\vec{y} = 1
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
This will generate two folders: one for the checkpoints during the finding of optimizations, and the other for the record of solution (might fail due to some reasons, such as memory limit or time limit, etc.). Now we should focus on the file `y.text` which records the components of the vector $\vec{y}$ in the SDP problem. From this file, we can reconstruct the original positivity bound $\vec{f}(x)\cdot\vec{y} \ge 0$ and check whether it is negative between some sampling points.

### Example
#### Problem description
Consider the positivity bound
$$
    1 + x^4 + y * ( x^4/12 + x^2 ) >= 0.
$$
Our goal is to maximize the quantity
$$
    -y.
$$
This is obviously a SDP problem. In terms of objection and normalization vector, we have
$$
    \vec{y} = (1, y), \quad \vec{f}(x) = ( f1[x], f2[x] ) = ( 1 + x^4, x^4/12 + x^2 ),
$$
$$
    \vec{b} = ( 0, -1 ) \text{, such that } \vec{b} \cdot \vec{y} = -y,
$$
and 
$$
    \vec{n} = ( 1, 0 ) \text{, such that } \vec{n} \cdot \vec{y} = 1.
$$
#### Sampling points
The degree of polynomials is $4$, so we choose $5$ sampling points. Note that we **treat them as constant functions at every sampling points**.
```python
    x = [0.06, 0.57, 1.63, 3.42, 6.49]
```
This generates $5$ SDP blocks (**constant constraint**):
$$
    f1[x1]*1 + f2[x1]*y \ge 0,
$$
$$
    f1[x2]*1 + f2[x2]*y \ge 0,
$$
$$
    f1[x3]*1 + f2[x3]*y \ge 0,
$$
$$
    f1[x4]*1 + f2[x4]*y \ge 0,
$$
$$
    f1[x5]*1 + f2[x5]*y \ge 0.
$$
#### SDP solution
These positive bounds give the solution
$$
    \vec{y} = ( 1, -2.483).
$$

HOWEVER, ploting the functional $\vec{f}(x) \cdot \vec{y}$, we find that it is negative around $0.75$ and $1.725$.

![image](/home/htk/sdpb_test/test1/different_y.png)

#### Manual Adaptive Sampling
Therefore, I take an adaptive sampling **manually**:
$$
    x = [0.06, 0.57, 1.0, 1.25, 1.50].
$$
and the corresponding solution is $\vec{y} = (1, -1.846)$, which is very close to the exact solution $(1, -1.840)$.

We see that adding the sampling points in the negative region of $f1(x) - 2.483 f2(x)$, the SDPB solver improves the solution a lot.

## Conculsion
**Adaptive sampling is necessary.** HOWEVER, there are some problems that cause this procedure to be time-consuming

1. We have to run the SDPB solver MANUALLY every re-sampling. 
2. After every SDPB solving, we have to MANUALLY check the solution in `y.txt`.
3. After checking the solution in `y.txt`, we need to verify the regions where the functional is negative, and consider the resampling again at these regions.

The STOPPING criterion is that the area (length) of the regions where the functional is negative is smaller then a given threshold. For example, size is $10^{-6}$ and the sampling is $x \in [x^* - Ns, x^* - (N-1)s, \cdots, x^*, \cdots, x^* + (N-1)s, x^* + Ns]$, where $x^*$ is the median of the negative region $[x_a, x_b]$ and $s = |x_b-x_a|/N$ and $N$ is set to $10$.

IMPROTANT: We have to iteratively write the PMP problems in to `JSON` files at different sampling points, and then turn them into `SDP` problems and then solve them by the SDPB solver. After each solving, we have to check the negative region of the functional with the solution, and determine the sampling points for the next run. 

It should be done by a `shell script` that has `sudo` permission.