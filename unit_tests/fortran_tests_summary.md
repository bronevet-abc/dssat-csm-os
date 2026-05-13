# Summary of Unit Tests for Biochar Fortran Module

This document summarizes the unit tests implemented for the Biochar module (`Biochar_mod.f90`) based on the hypotheses described in `hypothesis_tests.md`.

## Implemented Tests and Results

The following tests were implemented in [test_Biochar_mod.f90](file:///usr/local/google/home/kejie/Desktop/somic/dssat-csm-os/Soil/Biochar/test_Biochar_mod.f90) and successfully passed:

| Test # | Hypothesis / Name | Description | Result |
| :--- | :--- | :--- | :--- |
| **1** | The Nitrogen Cascade | Verified that an influx of ammonium accelerates the weathering of the ash pool (`M_BC`). | **PASSED** |
| **4** | Anion Competition | Verified that increasing Phosphate increases decomposition rate by reducing Nitrate sorption and increasing available Nitrogen. | **PASSED** |
| **5** | Dual-Driven Oxidation | Verified that oxidation rate decreases as potential CEC approaches its maximum. | **PASSED** |
| **7** | Precipitation vs Sorption | Verified that higher weathering (more Ca) increases P precipitation. | **PASSED** |
| **10** | Overflow Respiration | Verified that decomposition stops when Nitrogen pools are zero. | **PASSED** |
| **12** | Ash Exhaustion | Verified that weathering stops when the ash pool is exhausted. | **PASSED** |
| **14** | Feedback Cascade | Verified that high Labile Biochar suppresses ash weathering over a 24-hour integration due to N immobilization. | **PASSED** |
| **15** | Pure Thermal Oxidation | Verified that oxidation proceeds even without biological activity (recalcitrant C pool decay). | **PASSED** |
| **17** | Transient Immobilization Cliff | Verified that microbial immobilization collapses when the labile biochar pool is exhausted. | **PASSED** |

## Unimplemented / Skipped Tests

Some tests could not be fully implemented or were added as placeholders that report "SKIPPED" for the following reasons:

| Test # | Hypothesis / Name | Reason for Skipping |
| :--- | :--- | :--- |
| **2, 8** | Aluminium Thresholds | Require Aluminum penalty logic which is not present in the current `Biochar_mod.f90`. |
| **3, 19** | npSOM Bridging / AEC Erosion | Marked as placeholders or pending implementation in the source document. |
| **6, 18** | Plant Demand Competition | Require a plant N demand sink in the derivatives, which is noted as removed in the source document. |
| **9, 11** | Moisture/Temp Forcings | Require moisture or temperature effects on parameters (like pKa) that are currently static constants in the code. |
| **13** | Langmuir Exhaustion | Hard to verify without exposing internal state variables that are not returned by derivatives or updated in state. |
| **16** | Drought Freeze | Requires moisture dependence to arrest fluxes, which is missing. |
| **20** | Sink Overdraw | Requires plant N demand sink to test overdraw limits. |

## How to Run the Tests

To run the tests again, navigate to the directory and execute the test binary:

```bash
cd /usr/local/google/home/kejie/Desktop/somic/dssat-csm-os/Soil/Biochar
./test_Biochar_mod
```

If you need to recompile the tests, you can use the following command (assuming `.o` files are present or you compile them first):

```bash
gfortran -o test_Biochar_mod test_Biochar_mod.f90 CSMVersion.o OSDefsLINUX.o ModuleDefs.o Biochar_mod.o
```
