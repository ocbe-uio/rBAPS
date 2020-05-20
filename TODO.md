# Roadmap for next stable release

## Minimum requirements

For the first stable release of rBAPS, the following features should be implemented:

- [ ] Clustering of individuals (import `greedyMix.m` from BAPS)
- [ ] Clustering of populations (import `greedyPoPMix.m` from BAPS)
- [ ] Spatial clustering of individuals (import `spatialMixture.m` from BAPS)
- [ ] Spatial clustering of populations (import `spatialPoPMixture.m` from BAPS)
- [x] Admixture analysis (import `admix1.m` from BAPS)

Note to contributors: as tasks get finished (translated, not necessarily tested), please update this file with an `x` and keep it there. This should help us fill out a changelog for an eventual stable release.

## Wishlist

The list below contains non-essential but nice-to-have tasks for the next stable release.

- [ ] Implement sparse matrix optimizations from [`fastbaps`](https://github.com/gtonkinhill/fastbaps)
- [ ] Implement plotting functionality from [`starmie`](https://github.com/sa-lee/starmie)
- [ ] Standardize capitalization of function names
- [ ] Standardize capitalization of function arguments
- [ ] Replace redundant functions (ex.: `randga`)

# Known pitfalls

The following behavioral differences have been detected between the Matlab functions and their R counterparts. In order to save time, these differences will not be addressed, since they could require extensive reworking of a function. However, such differences may very well cause unexpected problems in some situations, which is why compiling this list is so important. The tables below might provide a good starting point for identifying and fixing bugs.

As general remarks, one should keep in mind that:

- For compliance with IEC 60559, the `round` in base R rounds .5 to the nearest even integer, whereas the homonym function in Matlab rounds up (or down, if negative).
- Some clobal variables have been added as a new (last) argument to the function they appear in.
- When a function is defined on Matlab with multiple outputs, as in `[y1, y2] = f(x)`, it will output only `y1` if called without assignment, as in `f(3)`, or if called with a one-length assignment, as in `a = f(3)`. In order to receive the full output, one must assign two variables to the left side of the assignment operator, as in `[a, b] = f(3)`. rBAPS Functions that might be affected by this behavior include `etsiParas`.

## Functional differences in rBAPS functions

| Function          | Argument           | Value          | Matlab output         | R output                        |
| ----------------- | ------------------ | -------------- | --------------------- | ------------------------------- |
| `ownNum2Str`      | `number`           | `NaN`          | `'NAN'`               | error                           |
| `ownNum2Str`      | `number`           | `<vector>`     | `'<vector elements>'` | `'<vector elements>'` + warning |
| `repmat`          | `length(n)`        | `> 2`          | > 2D matrix           | 2D matrix                       |

## Functional differences in base Matlab functions

Function | Matlab output | R output
-------- | ------------- | --------
`max` | Can handle complex numbers | Cannot handle complex numbers