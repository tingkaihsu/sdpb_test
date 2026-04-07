# Discretization (Adaptive Refinement)
SDP problem
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
