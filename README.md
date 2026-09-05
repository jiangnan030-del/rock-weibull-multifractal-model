# Rock strength Weibull model from multi-fractal damage

This repository is a Python implementation of the paper currently open in Notion:

**刘树新 等 (2011)** — *基于损伤多重分形特征的岩石强度 Weibull 参数研究*.

## What is implemented

The code converts the core equations in the paper into reusable Python functions:

- micro-element strength based on the Mohr-Coulomb criterion (Eq. 8)
- Weibull damage variable (Eq. 7)
- constitutive relation with fractal damage parameter `lambda` (Eq. 9)
- closed-form solution of Weibull parameters `m` and `F0` from the peak-point conditions (Eqs. 18-19)
- stress-strain curve generation for different crack-spacing / multifractal cases
- reproduction datasets for Table 3 and Table 4
- example figures generated with matplotlib

## Files

- `rock_weibull_model.py`: main implementation
- `data/table3_relative_dimension.csv`: digitized Table 3 data
- `data/table4_singularity.csv`: digitized Table 4 data
- `figures/`: generated example plots

## Notes

The paper page provides the governing formulas and summary tables, but not a complete experimental parameter set for every figure. Therefore this repo includes:

1. exact implementations of the equations;
2. digitized values from Tables 3 and 4;
3. a configurable demo parameter set for example stress-strain curves.

If you later provide a specific rock's experimental parameters (`E`, `nu`, `phi`, `sigma1c`, `epsilon1c`, `sigma3`), you can plug them directly into `PeakState` and regenerate case-specific results.

## Quick start

```bash
python3 rock_weibull_model.py
```

This writes CSV tables and PNG figures into the repository.

## Example usage

```python
from rock_weibull_model import PeakState, solve_weibull_parameters

peak = PeakState(
    E=15000.0,
    nu=0.25,
    sigma3=0.0,
    sigma1c=35.0,
    epsilon1c=0.004,
    phi_deg=35.0,
)

params = solve_weibull_parameters(peak, lambda_value=0.9293)
print(params)
```

## Dependencies

- Python 3.10+
- matplotlib

## Possible next steps

- fit model parameters against your own stress-strain test data;
- add a Jupyter notebook for parameter inversion;
- implement direct multifractal-spectrum calculation from crack-network data;
- compare relative-dimension and singularity-based `lambda` choices.
