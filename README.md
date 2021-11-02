
<!-- README.md is generated from README.Rmd. Please edit that file -->

# relapseRisk

<!-- badges: start -->
<!-- badges: end -->

R package for predicting risk of relapse for anorexia nervosa recovery.
UNDER CONSTRUCTION.

## Background

The following documents demonstrate analyses involved in the
construction and evaluation of the machine learning algorithm used to
predict recovery outcomes in patients diagnosed with anorexia nervosa.
To view these files, the user must download the raw files and then open
them in a local web browser.

[EPSI\_thetaAnalysis](./EPSI_thetaAnalysis.html) contains an overview of
the initial analysis that led to the construction of this package. This
involved analyzing data from previous research on eating disorder
patients who completed the Eating Pathology Symptoms Inventory (EPSI)
during treatment, and were then re-evalued by therapists during a
6-month and 12-month follow-up to determine whether they had relapsed
based on their body-mass index (BMI), number of binge-eating episodes
(BE), and compensatory episodes (CE). The objective was to determine
whether higher-order factors derived from the patients’ EPSI responses
were useful in predicting their recovery outcomes.

[fullAnalysis](./fullAnalysls.html) contains an analysis of the training
data used to construct the machine learning algorithm that predicts the
risk of relapse for new patients diagnosed with anorexia nervosa. The
purpose of this analysis was to evaluate multiple resampling algorithms
and performance outcomes to select the model that best predicts recovery
outcomes within the training data. The model selected based on these
results will then be used in our current research to predict patients’
risk of relapse over the course of a 12-week randomized controlled trial
(RCT). Upon re-assessing patients at a 6-month follow-up, results are
then aggregated with the training data to further refine the algorithm
and improve prediction for future patients.
