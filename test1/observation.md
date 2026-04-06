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


## 3rd Update
### Previous Inconsistency
The auto-built version is in fact still analytic form: 

```
[ [ [ [ "1", "0", "0", "0", "1" ], [ "0", "0", "1", "0", "0.0833"] ] ] ]
```
that correspond to $f1[x] = 1 + x^4$ and $f2[x] = x^2 + x^4/12$. The `polynomial` IN FACT stores the coefficient of $f1[x]$ and $f2[x]$.

### Corrections

For non-polynomial functions, I initially plan to **Talyor Expand** the functions and extract the coefficients of each term.

My collaborator suggests a different approach: treat them as CONSTANT.

Now PMP per sample point, which is weird at first. 

The structure is now:
A PMP for a sample point
```
{
    "samplePoints":[
        "x1"
    ],
    "sampleScalings":[
        "s1"
    ],
    "polynomials":[
        [
            [
                [
                    [
                        "f1[x1]"
                    ],
                    [
                        "f2[x1]"
                    ]
                ]
            ]
        ]
    ]
}
```
for next sample point, we open another pmp
```
{
    "samplePoints":[
        "x2"
    ],
    "sampleScalings":[
        "s2"
    ],
    "polynomials":[
        [
            [
                [
                    [
                        "f1[x2]"
                    ],
                    [
                        "f2[x2]"
                    ]
                ]
            ]
        ]
    ]
}
```
etc.


### Result
My result: -2.483

The analytic result: -1.840

My result with different sample point: -1.846
The analytic result: -1.840
