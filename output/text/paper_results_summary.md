# Summary of Key Findings

**Headline:** The highest healing probability was observed for acropora treated with dremel (35.1%; 95% CI [0.21, 0.52]).


## Model comparison
- Treatment × species interaction: χ²(2) = 3.553, p = 0.169.
  (From `anova(mod_add, mod_int)`.)
## Main effects (drop1 on mod_add)
- treatment: χ²(1) = 16.248, p < 1e-4
- species: χ²(2) =  4.803, p = 0.091

## Marginal probabilities (emmeans; response scale)
- Highest Pr(healed): acropora × dremel: 35.1% (95% CI [0.21, 0.52])
- Lowest  Pr(healed): pocillopora × airbrush: 0.4% (95% CI [0, 0.04])
## Treatment contrast within species (Dremel vs Airbrush)
- Within acropora, Dremel had higher Pr(healed) vs Airbrush by 33.1% (95% CI 17.3%, 48.8%); OR = NA (95% CI NA, NA), p < 1e-4
- Within porites, Dremel had higher Pr(healed) vs Airbrush by 19% (95% CI 2.2%, 35.9%); OR = NA (95% CI NA, NA), p = 0.027
- Within pocillopora, Dremel had higher Pr(healed) vs Airbrush by 9.2% (95% CI 2.9%, 21.3%); OR = NA (95% CI NA, NA), p = 0.138

## Species contrasts within treatment (Tukey-adjusted)
- Under dremel, acropora vs pocillopora in Pr(healed) was 25.5% (95% CI 1.7%, 49.4%); OR = 5.11 (95% CI 0.75, 34.65), p = 0.032
- Under dremel, acropora vs porites in Pr(healed) was 15.1% (95% CI 12.6%, 42.8%); OR = 2.16 (95% CI 0.47, 9.96), p = 0.408
- Under airbrush, acropora vs pocillopora in Pr(healed) was 1.6% (95% CI 2.3%, 5.6%); OR = 5.11 (95% CI 0.75, 34.65), p = 0.599
- Under dremel, pocillopora vs porites in Pr(healed) was 10.4% (95% CI 36.2%, 15.3%); OR = 0.42 (95% CI 0.05, 3.72), p = 0.609
- Under airbrush, acropora vs porites in Pr(healed) was 1.1% (95% CI 2%, 4.2%); OR = 2.16 (95% CI 0.47, 9.96), p = 0.687
- Under airbrush, pocillopora vs porites in Pr(healed) was 0.5% (95% CI 2.4%, 1.3%); OR = 0.42 (95% CI 0.05, 3.72), p = 0.773

## Model diagnostics
- Additive model (mod_add): overdispersion = 1.107, singular = TRUE, R²_m = NA, R²_c = NA, AIC = 99.35, BIC = 117.45.
- Interaction model (mod_int): overdispersion = 0.815, singular = TRUE, R²_m = NA, R²_c = NA, AIC = 99.80, BIC = 123.94.
- Latent-scale ICCs: parent_id = 0.000, coral_id-within-parent = 0.000.
