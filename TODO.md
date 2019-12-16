# Minimum requirements for next stable release

For the first stable release of rBAPS, the following features should be implemented:

- [ ] Clustering of populations (import `greedyPoPMix.m` from BAPS)
- [ ] Clustering of individuals (import `greedyMix.m` from BAPS)
- [ ] Spatial clustering of populations (import `spatialPoPMixture.m` from BAPS)
- [ ] Spatial clustering of individuals (import `spatialMixture.m` from BAPS)
- [ ] Admixture analysis (import `admix1.m` from BAPS)

Note to contributors: as tasks get finished, please update this file with an `x` and keep it there. This should help us fill out a changelog for an eventual stable release.

# Known pitfalls

The following behavioral differences have been detected between the Matlab functions and their R counterparts. In order to save time, these differences will not be addressed, since they could require extensive reworking of a function. However, such differences may very well cause unexpected problems in some situations, which is why compiling this list is so important. The list below might provide a good starting point for identifying and fixing bugs:

## `ownNum2Str`

Argument | Value | Matlab output | R output
---------|-------|---------------|---------
`number` | `'NaN` | `'NAN'` | error
`number` | `<vector>` | `'<vector elements>'` | `'<vector elements>'` + warning