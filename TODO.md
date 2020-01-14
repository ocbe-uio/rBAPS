# Roadmap for next stable release

## Minimum requirements

For the first stable release of rBAPS, the following features should be implemented:

- [ ] Clustering of populations (import `greedyPoPMix.m` from BAPS)
- [ ] Clustering of individuals (import `greedyMix.m` from BAPS)
- [ ] Spatial clustering of populations (import `spatialPoPMixture.m` from BAPS)
- [ ] Spatial clustering of individuals (import `spatialMixture.m` from BAPS)
- [ ] Admixture analysis (import `admix1.m` from BAPS)

Note to contributors: as tasks get finished, please update this file with an `x` and keep it there. This should help us fill out a changelog for an eventual stable release.

## Wishlist

The list below contains non-essential but nice-to-have tasks for the next stable release.

- [ ] Uniformize capitalization of function names
- [ ] Uniformize capitalization of function arguments
- [ ] Implement optimizations from `fastbaps`
- [ ] Replace redundant functions (ex.: `randga`)

# Known pitfalls

The following behavioral differences have been detected between the Matlab functions and their R counterparts. In order to save time, these differences will not be addressed, since they could require extensive reworking of a function. However, such differences may very well cause unexpected problems in some situations, which is why compiling this list is so important. The list below might provide a good starting point for identifying and fixing bugs:

Function | Argument | Value | Matlab output | R output
---------|----------|-------|---------------|---------
`ownNum2Str` | `number` | `NaN` | `'NAN'` | error
`ownNum2Str` | `number` | `<vector>` | `'<vector elements>'` | `'<vector elements>'` + warning
`repmat` | `length(n)` | `> 2` | > 2D matrix | 2D matrix

As general remarks, one should keep in mind that:

- For compliance with IEC 60559, the `round` in base R rounds .5 to the nearest even integer, whereas the homonym function in Matlab rounds up (or down, if negative).
- Some clobal variables have been added as a new (last) argument to the function they appear in.