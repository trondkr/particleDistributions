# particleDistributions
A toolbox for creating distributional figures from particle tracking simulations using OpenDrift.

The following calculations needs to be added:

- [ ] Sedimentation maps that correctly calculates the sedimentation timestamp based on the distance from bottom depth and the average movement of the particle over the last 2 hours. Sedimentation maps needs to be created for each week of the year combining all polygons. Should include all particle types that have sedimented.

- [ ] Create TIFF files from the sedimentation maps (for import to ARCGIS)

- [ ] Calculate statistics for each polygon - a histogram showing percentages of particles that have sedimented for each particle type and the distance sedimented from polygon.

- [ ] Finish the vertical profile plot for individual particle sizes combining all polygons, and removing the data after particle has sedimented.

- [ ] The sedimentation function should take distance from bottom and maximum movement over X hours as input to determine when sedimentation occurss to make it more generic.

- [ ] Calculate the probability that kelp from a given polygon is sedimented at depth deeper than X meters.


