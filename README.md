# AutoEmma
Automated Efficient Mixed Model Association Studies - R Package

## Installation
1. Download AutoEmma_1.2.1.tar.gz from GitHub.
2. Open R-Studio and type the following command into the Console.
```
install.packages("file/path/AutoEmma_1.2.1.tar.gz", repos = NULL, type = "Source")
```

## Quick Tutorial
1. Load AutoEmma Package into your current R session: ```library(AutoEmma)```
2. Create diagnostic graphs for your input files: ```ae.input.graphs(gfile = "your_genotypes.csv", pfile = "your_phenotypes.csv")```
3. Create a map file: ```my_map <- ae.make.map(file = "your_phenotypes.csv")```
4. Run AutoEmma: ```results <- ae.run(genotype = "your_genotypes.csv", phenotype = "your_phenotypes.csv", map = my_map, method = "REML")```
5. Customize graphs if desired: ```ae.manhattan(df = results, file_name = "new_manhattan.png", point_size = 1.2, colors = c("royalblue3", "gray18"), title = "Manhattan Plot", sl = 4.5, pval = 0.05, snps = 100000, exclude = 0.5, annotate = TRUE)```

## Change Log
**Version 1.0.0**

* Initial release

**Version 1.0.1**

* Simplified input diagnostic graphs
* Fixed a bug in strain selection for cladogram creation

**Version 1.1.1**

* Improved manhattan and Q-Q plots
* Added annotation of suggestive SNPs to manhattan plot

**Version 1.2.0**

* Added SEM option to average phenotype plot
* Improved look of input diagnostics graphs
* Vectorized averages calculation to speed up input diagnostics graphs

**Version 1.2.1**

* Added option to choose between wide-screen graphs (16x9) and square graphs (9x9)
* Improved text updates for input diagnostic graphs
