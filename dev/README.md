
## Generating the `mediation_example` data

Initializing data were downloaded from the
[zip](https://cdn1.sph.harvard.edu/wp-content/uploads/sites/603/2019/03/MediationPsychMethods.zip)
available at Tyler VanderWeele’s
[site](https://www.hsph.harvard.edu/tyler-vanderweele/tools-and-tutorials/)
(accessed 2020-04-13). This data set has 100 observations. We resampled
with replacement to boost the sample size to 250 to develop and test a
variety of features in the `mediator` package.

``` r
# load VanderWeele data
mediation_survival <- readr::read_delim(
  here::here("dev/example250/mediation_survival.txt"),
  delim = " ",
  skip = 4,
  n_max = 100,
  col_names = c("id", "x", "m", "y", "cens", "c")
)

# bind 150 randomly resampled rows
set.seed(8675309)
mediation_example <- 
  dplyr::bind_rows(
    mediation_survival,
    dplyr::sample_n(
      mediation_survival,
      150,
      replace = TRUE)
    )
  
# select and rename variables
mediation_example <- 
  mediation_example %>% 
  dplyr::select(x, m, y, cens, c) %>%
  dplyr::rename(
    m_01 = m,
    m = y,
    y = cens
    )
    
# store
readr::write_rds(
  mediation_example, 
  here::here("data/mediation_example.rds")
)  
```
