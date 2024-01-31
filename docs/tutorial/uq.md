<!-- 
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is under the MIT license found in the root directory LICENSE.md 
Everyone is permitted to copy and distribute verbatim copies of this license 
document.

This version of the MIT Public License incorporates the terms
and conditions of MIT General Public License.
-->
## Uncertainties in Underwater GEFs and Induced Submarine Cable Induced Voltages
The primary goal is to employ advanced techniques, including Monte Carlo simulations, to comprehensively assess the uncertainties associated with GEFs and induced voltages in submarine cables. By providing a nuanced understanding of these uncertainties, the study aims to contribute to the development of more accurate and reliable models for assessing the impact of Geomagnetic Disturbances on critical underwater infrastructure. The following paragraphs shows step-by-step process that leads to `UQ` in GMD-induced GEFs and cable voltages. Furthermore, this tool is also used to conduct **Sensitivity Analysis**.

### Understanding Cable Route Dynamics:
> When calculating induced voltages along submarine cable routes, it is crucial to account for variations in Ocean-Earth thickness structures and the number of cable sections. Existing studies, including Boteler et al. (2003; 2019) and Chakraborty et al. (2022), emphasize the significant impact of cable section length and water depth on calculated voltages.

### Consideration of Uncertainties:
> In addition to these structural variations, uncertainties in magnetic field values used for voltage calculations are paramount. Interpolated geomagnetic field values can introduce uncertainties, making it essential to employ methods that consider these variations.

### Addressing Structural Variations:
> To tackle variations in the Ocean-Earth thickness structure within cable sections and uncertainties in interpolated magnetic field values, a robust approach involves characterizing these variations. Commonly used techniques include sensitivity analysis, Monte Carlo simulations, and uncertainty propagation methods.

### Enhancing Reliability Through Uncertainty Analysis:
> Incorporating uncertainty analysis significantly improves the reliability of voltage calculations and model predictions. By acknowledging and quantifying uncertainties, the study aims to provide a more comprehensive understanding of variations in Ocean-Earth thickness structures and induced voltage outcomes.

### Monte Carlo Simulations:
> The study will leverage Monte Carlo simulations, specifically Markov Chain Monte Carlo (MCMC), to quantify variations in Ocean-Earth thickness structures and uncertainties in interpolated magnetic field estimates. This method involves sampling input parameters from distributions and propagating them through the SCUBAS model. The specific steps involves:
>
> * Divide the cables into different segments $$\mathcal{N}$$ (e.g., 9).
> * For each segment create virtual 10-(user provided $$\mathcal{V}$$) samples of latitude/longitude points along the cable (say, segment-sample-virtual-point `SSVP`).
> * For $$\mathcal{N}$$ and $$\mathcal{V}$$ `SSVP` we have $$\mathcal{N}\times\mathcal{V}$$ possible Ocean-Earth structures to run simulations. 
> * Run SCUBAS for all these $$\mathcal{N}\times\mathcal{V}$$ combinations, and from these solutions quantify uncertainty in the output use MC-simulation.
> * For now we have no-way to estimate error in input magnetic field dataset. If given we can generate another set of magnetic field input and then run MC-simulation to propagate the errors in the output.


### Statistical Representation of System Response:
> Monte Carlo simulations generate numerous possible outcomes, offering a statistical representation of the system's response under varying conditions. Uncertainties in the Ocean-Earth thickness structure will be addressed by repeatedly sampling ocean water levels from each cable section along the route.

### Characterizing Uncertainty:
> Each simulation run produces a distinct set of cable voltages, resulting in an output response distribution. Statistical measures such as mean, standard deviation, confidence intervals, or probability density functions will characterize the uncertainty in GMD-induced cable voltage.

### Insights and Sensitivity Analysis:
> Monte Carlo simulations not only provide insights into the range of system behavior but also enable sensitivity analysis. This analysis allows assessment of the influence of different input parameters on output variability, contributing to a more nuanced understanding of induced voltage outcomes under various conditions.