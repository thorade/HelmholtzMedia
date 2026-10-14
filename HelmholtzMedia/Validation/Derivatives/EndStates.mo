within HelmholtzMedia.Validation.Derivatives;
type EndStates = enumeration(
    SinglePhase "supercritical states",
    TwoPhase "two-phase states",
    PhaseChange "from subcooled liquid through the two-phase region to superheated vapour",
    LiquidNearSaturation "liquid states, ending close to the saturated liquid",
    VapourNearSaturation "vapour states, starting close to the saturated vapour")
  "end states of the linear ramps in CheckDerivativeAnnotations";
