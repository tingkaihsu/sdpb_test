# Observation
1. New reducedPrefactor in the numeric implementation
2. Polynomials with wrong format in the numeric implementation
3. Missing reducedSampleScalings, bilinearBasis_0, and bilinearBasis_1
   in the new numeric implementation.
## Update 
1. New reducedPrefactor in the numeric_pmp.json (minor)
2. bilinearBasis_0, reducedSampleScalings, and bilinearBasis_1 STILL missing (minor)
3. polynomials format is STILL wrong. (**Important**)

Polynomials format in numeric implementation:
[ [ f1[x1], f2[x1] ], [ f1[x2], f2[x2] ], ..., [ f1[x5], f2[x5] ] ].

The correct format:
[ [ [ [ f1[x1], ..., f1[x5] ], [ f2[x1], ..., f2[x5] ] ] ] ].


