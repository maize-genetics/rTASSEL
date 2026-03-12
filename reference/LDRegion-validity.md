# LDRegion validation

Validates that an `LDRegion` object has properly formed slots: `start`
and `end` must be length-1 finite numerics with `end >= start`, `color`
must be a recognized R color, `linewidth` must be length-1 (positive or
NA), and `showSpan` must be a length-1 non-NA logical.

## Arguments

- object:

  An `LDRegion` object.
