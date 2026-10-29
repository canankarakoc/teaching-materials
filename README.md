# Teaching materials

Hands-on activities, worksheets and R scripts I designed for courses in ecology, evolution and mathematical biology. Each folder is one course; every script runs on its own from its folder.

## Courses

### Biodiversity lab course, UFZ Leipzig (2018) · `biodiversity-lab-2018/`

Students sample gummy-candy "communities" of different sizes, then calculate diversity indices on their counts and on 18S fingerprinting (T-RFLP) data from the course.

- `diversity_exercise.R`: richness, Shannon and Simpson diversity and their Hill-number forms, relative abundances and NMDS ordination (`vegan`)
- `gummy_species_key.pdf`: species key for the gummy candies
- `HariboMatrix.csv`: gummy counts per community and sample size
- `ClusBinMatrix_B_3_Praktikum_18S.txt`: 18S T-RFLP fragment matrix

### Practical training: traits and fitness costs of evolved populations, Leipzig University (2020) · `experimental-evolution-2020/`

BSc Biology, module 11-BIO-0630 (Environmental Microbiology), run at UFZ. Students compete evolved *E. coli* clones against their ancestor and measure growth and phage infectivity, to estimate the fitness costs of resistance.

- `Practical_training_protocol_2020` (PDF and Word): course protocol with the two-week schedule, background, methods and media recipes
- `Measuring_fitness_costs` (PDF and Word): background and methods
- `competition/class.R`: in-class walk-through of the competition data
- `competition/fitness_costs.R`: selection coefficients and relative fitness
- `growth/class_14.05.R`: growth-curve fitting (`growthrates`)
- `competition_data.csv`, `growth_data2.csv`, `plate_design2.csv`: course data

### Quantitative Biodiversity (Z620), Indiana University (2023) · `quantitative-biodiversity-2023/`

My contributions to the course: sample coverage and Hill numbers, and a gummy-candy sampling exercise with a matching simulation. Each `.Rmd` comes with its rendered PDF.

- `Introduction_iNEXT`: equalizing samples by coverage, Hill numbers and extrapolation (`iNEXT`)
- `GummyCandyExercise`: three gummy species pools with different richness and evenness, sampled by groups
- `Simulations_of_sampling`: simulated communities and sampling for the same exercise (`mobsim`, `vegan`, `iNEXT`)

### BIOL 3400, Georgia Institute of Technology (2025) · `mathbio-2025/`

- `Probability2_Bayes_worksheet` (PDF and Word): Bayes' theorem with diagnostic tests (HIV testing, sequential updating, mammography)
- `Beans_coins_dice_activity` (PDF and Word): students build a 2×2 diagnostic table by hand, using beans for disease status and a coin and a die as the test
- `sequential_updating.nb`: sequential Bayesian updating (Mathematica)
- `Selection_classExercise.nb`: allele frequency change under haploid selection (Mathematica)
- `SIR_simulation.R`: SIR epidemic model (`deSolve`)

## How to use

Set the working directory to the folder a script is in, then run it top to bottom. Each script reads only files from its own folder.

Last run in October 2026 with R 4.6.0, vegan 2.7.6, ggplot2 4.0.3, reshape 0.8.10, tidyverse 2.0.0, growthrates 0.8.5, deSolve 1.42, scatterplot3d 0.3.45, iNEXT 3.0.2 and mobsim 0.3.2.

```r
install.packages(c("vegan", "ggplot2", "reshape", "tidyverse", "growthrates", "deSolve", "scatterplot3d", "iNEXT", "mobsim"))
```

The `.nb` notebooks open in Mathematica or in the free Wolfram Player.

Answer keys are not included. Instructors can ask me for them.

## License

Code (`.R`, the code in `.Rmd` files, and `.nb` notebooks) is under the MIT License (`LICENSE`). Worksheets, handouts, images and data are under CC BY 4.0 (`LICENSE-CC-BY-4.0.md`).
