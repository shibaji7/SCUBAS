# Component: transfer function

:::{note}
Draft stub. The writer filled in what the code and the papers say; the authors complete the physics below. Section 3 ("To fill in") is yours.
:::

## 1. In the code today

- `scubas.models.OceanModel`: `calcZ` (impedance of the layered Earth for a grid of frequencies, default 1e-6 to 1 Hz), `calcTF` and `get_TFs` (the transfer functions `E2B` and `B2B`),
  `to_Efields` (FFT of the record, multiplication, inverse FFT). `Preprocess` tapers and detrends the record before the transform.
- The definitions of `E2B` and `B2B` as the code writes them are in the [glossary](../../glossary.md).
- The legacy constants (taper parameter, baseline length, frequency limits) are listed with their status in the [parameter list](../../dev/parameters.md).

## 2. In the papers

- The transfer function relates the seafloor electric field to the surface magnetic field and depends on the seawater depth and conductivity and on the surface impedance of the Earth below
  (Boteler et al., 2024, Section 3.1, Equations 13 and 14, from Boteler and Pirjola, 2003). Its amplitude falls with water depth, and its phase changes sign with frequency for deep
  water (Chakraborty et al., 2026, Section 2.3).
- The frequency-domain relations are $E_N = T_X B_E$ and $-E_E = T_X B_N$ (Equation 16).
- The skin depth in the Earth is less than a sixth of a wavelength, so the induced current is concentrated near the surface and there is no return current at depth (Boteler et al., 2024, Section 2.1).

## 3. To fill in (authors)

- The full derivation of the transfer function in the form the code uses, and how it maps to Equation 13 of Boteler et al. (2024).
- The frequency range and taper, and the reason for them.
- The effect of the record's sampling and length on the result.

Back to [the model](../index.md).
