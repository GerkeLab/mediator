
# mediator <img src="man/figures/hex.png" align="right" height="139" />

The goal of mediator is to conduct causal mediation analysis under the
counterfactual framework, allowing interation between the exposure and
mediator (cite VanderWeele here). (State the kinds of effects that are estimated: NDE, CDE, proportion mediated, etc).

## Installation

You can install `mediator` from github with:

``` r
# install.packages("devtools")
devtools::install_github("gerkelab/mediator")
```

## Usage

`mediator` currently implements mediation analysis for dichotomous and
count mediators and outcomes (also censored time-to-event outcomes now?). (I think a table would be helpful here to summarize the different kinds of variables and model types that `mediator` fits). Estimate validity assumes proper modeling on the part of
the user.

(Suggestion: I would place a DAG here for some type of mediation mechanism, then state "we want to estimate the direct and indirect effect of _whateverVariable_ on _whateverOutcomeVariable_". You can then demo copy-paste-able code segments with output and an explanation of how we interpret that output with respect to our original question.)

## Additional resources

For an in-depth explanation of mediation analysis or complementary tools
for SAS or SPSS users, please check out Linda Valeri and Tyler
VanderWeele’s paper and macros, which are available on VanderWeele’s
[website](https://www.hsph.harvard.edu/tyler-vanderweele/tools-and-tutorials/).

The parametric model-based approach of `mediator` differs from another R package which conducts mediation analyses under a non-parametric framework, [`mediation`](https://cran.r-project.org/web/packages/mediation/vignettes/mediation.pdf).
