# Parameter Uncertainty in Epidemiological Modeling using Generalized Growth Model (GGM)

This project explores the estimation and uncertainty analysis of epidemiological parameters using the Generalized Growth Model (GGM). We fit the model on both synthetic and real epidemic data, and use bootstrapping to quantify parameter uncertainty and the effective reproduction number \( R_t \).

## 📁 Project Structure

. ├── data/ # Contains synthetic and real epidemic datasets ├── images/ # Visualizations (model fit, bootstrap histograms, etc.) ├── scripts/ # MATLAB code for model fitting, bootstrapping, and plotting ├── IP_Project_Parameter_Uncertainty_on_Epidemiology.pdf # Final report ├── README.md # This file └── .zenodo.json # (Optional) Metadata for Zenodo archiving


## 📌 Simulations Overview

1. **Simulation 1** – Fitting GGM to simulated exponential growth data  
2. **Simulation 2** – Fitting GGM to synthetic data generated using GRM  
3. **Simulation 3** – Bootstrap estimation of uncertainty for synthetic GRM data  
4. **Simulation 4** – Model fit and forecasts from bootstrap samples  
5. **Simulation 5** – Application of GGM to early-phase data from the 1918 Influenza Pandemic in San Francisco

Each simulation involves:
- Estimation of GGM parameters \( r \) (growth rate) and \( p \) (deceleration)
- Computation of time-varying effective reproduction number \( R_t \)
- Bootstrap-based uncertainty quantification
- Visualizations of model fits and parameter distributions

## 📊 Key Concepts

- **GGM (Generalized Growth Model):**
  \[
  \frac{dC(t)}{dt} = rC(t)^p
  \]
  Where:
  - \( C(t) \): cumulative incidence
  - \( r \): growth rate
  - \( p \): deceleration parameter

- **Bootstrap Realizations:**  
  200 resamples used to estimate variability in \( r, p, R_t \)

- **Effective Reproduction Number \( R_t \):**  
  Estimated using the incidence time series and a discretized exponential generation interval.

## 📄 Final Report

The complete methodology, figures, equations, and interpretation are documented in  
📌 [`IP_Project_Parameter_Uncertainty_on_Epidemiology.pdf`](./IP_Project_Parameter_Uncertainty_on_Epidemiology.pdf)

## 🛠 Tools & Environment

- **Language:** MATLAB  
- **Main Libraries:** `lsqcurvefit`, `bootstrp`, plotting utilities  
- **Visualization:** Histograms, line plots with confidence intervals

## 🧠 Authors

- Yash Goyal (Primary Contributor)

## 📜 License

This project is licensed under the MIT License.

## 📚 References

- Chowell et al., “A Novel Sub-Epidemic Modeling Framework for Short-Term Forecasting Epidemic Waves”  
- Viboud et al., “The RAPIDD ebola forecasting challenge: Synthesis and lessons learned”

---

