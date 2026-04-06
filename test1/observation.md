# Observation

## 2nd Update
1. Now the polynomials format is correct. 
2. But the problem is the values of polynomials are different, but they IN FACT correspond to same optimization problem.

For $f1[x] = 1+x^4$ and $f2[x] = x^4/12 + x^2$,
The original implementation:
```
[ [ [ [ "1", "0", "0", "0", "1" ], [ "0", "0", "1", "0", "0.0833"] ] ] ]
```
The numeric implementation:
```
[ [ [ [ "1.0", "1.1", "8.1", "140.", "1800." ], [ "0.0036 "0.33", "3.2", "23.", "190."] ] ] ]
```

Which are very different. I have checked the sampling points and the corresponding f1[x] and f2[x] values, which are correct. 

*Comment*: maybe some other modifications is used.
