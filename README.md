# structural-equation-modeling-latent-variables-analysis
Statistical modeling using structural equation models (SEM) for latent variable analysis with reliability, validity, and hypothesis testing.

# Structural Equation Modeling for Latent Variable Analysis

## Introduction

Understanding complex behavioral and organizational phenomena often requires modeling relationships between variables that are not directly observable. These constructs, known as latent variables, must be inferred from observed indicators.

From a statistical perspective, this leads to the use of **Structural Equation Modeling (SEM)**, a framework that combines multivariate analysis, regression, and factor analysis.

From a data science standpoint, SEM can be interpreted as a **structured predictive modeling framework**, where latent representations are learned and relationships between variables are quantified.

This project applies SEM to analyze relationships between latent constructs, combining:

- Statistical inference
- Measurement validation
- Predictive modeling concepts

---

## Dataset and Variables

The dataset consists of $n = 228$ observations, including:

- Observed indicators for latent constructs
- Sociodemographic variables
- Contextual variables

Exploratory data analysis and preprocessing were conducted prior to modeling :contentReference[oaicite:0]{index=0}.

---

## Statistical Methodology

### Nonparametric Hypothesis Testing

To compare group differences in median scores, nonparametric tests were applied:

- Mann–Whitney test (two groups)
- Kruskal–Wallis test (three or more groups)

The test statistic is based on rank sums:

$$
H = \frac{12}{N(N+1)} \sum_{g=1}^{G} n_g \bar{R}_g^2 - 3(N+1)
$$

where:

- $n_g$ is the sample size in group $g$
- $\bar{R}_g$ is the average rank in group $g$

Significant differences were identified for variables such as municipality size and demographic factors :contentReference[oaicite:1]{index=1}.

---

## Measurement Model

Latent variables are defined through observed indicators using a factor model:

$$
\mathbf{x} = \Lambda \boldsymbol{\xi} + \boldsymbol{\delta}
$$

where:

- $\mathbf{x}$: observed variables  
- $\boldsymbol{\xi}$: latent variables  
- $\Lambda$: factor loadings  
- $\boldsymbol{\delta}$: measurement error  

### Reliability and Validity

The following metrics were used:

#### Cronbach's Alpha

$$
\alpha = \frac{k}{k-1} \left(1 - \frac{\sum \sigma_i^2}{\sigma_T^2}\right)
$$

#### Composite Reliability

$$
CR = \frac{(\sum \lambda_i)^2}{(\sum \lambda_i)^2 + \sum \theta_i}
$$

#### Average Variance Extracted

$$
AVE = \frac{\sum \lambda_i^2}{\sum \lambda_i^2 + \sum \theta_i}
$$

Results indicated:

- $\alpha > 0.70$
- $CR > 0.70$
- $AVE > 0.50$

confirming internal consistency and convergent validity :contentReference[oaicite:2]{index=2}.

---

## Structural Model

The relationships between latent variables are modeled as:

$$
\boldsymbol{\eta} = B\boldsymbol{\eta} + \Gamma \boldsymbol{\xi} + \boldsymbol{\zeta}
$$

where:

- $\boldsymbol{\eta}$: endogenous latent variables  
- $\boldsymbol{\xi}$: exogenous latent variables  
- $B$: relationships among endogenous variables  
- $\Gamma$: effects of exogenous variables  
- $\boldsymbol{\zeta}$: structural errors  

---

## Moderation Effects

Interaction effects were included:

$$
\text{Intention} = \beta_1 A + \beta_2 N + \beta_3 PBC + \beta_4 (A \times PBC) + \beta_5 (N \times PBC)
$$

where:

- $A$: Attitude  
- $N$: Perceived Norm  
- $PBC$: Perceived Behavioral Control  

No significant moderation effects were observed :contentReference[oaicite:3]{index=3}.

---

## Model Evaluation

Model performance was assessed using:

### Coefficient of Determination

$$
R^2 = 1 - \frac{\text{Var}(\varepsilon)}{\text{Var}(Y)}
$$

Results:

- $R^2 \approx 0.54 - 0.72$ for intermediate constructs  
- $R^2 \approx 0.56 - 0.58$ for intention  

indicating good explanatory power :contentReference[oaicite:4]{index=4}.

---

## Results

Key findings:

- Significant relationships between latent variables  
- Strong measurement reliability and validity  
- Attitude significantly influences intention  
- No evidence of moderation effects  

From a data science perspective:

- The model explains up to ~58% of variance  
- Latent feature modeling improves interpretability  
- SEM acts as a structured predictive system  

---

## Conclusion

This work demonstrates how SEM provides a powerful framework for:

- Modeling latent variables  
- Testing theoretical relationships  
- Combining inference and prediction  

By integrating statistical rigor with data-driven modeling, this approach bridges traditional statistical analysis and modern data science workflows.

---

## References

- Chin, W. W. (1998). Partial Least Squares SEM  
- Cronbach, L. J. (1951). Coefficient Alpha  
- Fornell, C., & Larcker, D. (1981). Validity in SEM  
- Hair, J. F. et al. (2010). Multivariate Data Analysis  
