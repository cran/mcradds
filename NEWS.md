# mcradds 1.0.1

### Meta
* Remove `mcr` package related codes as it's not available in the CRAN.

# mcradds 1.0.0

### Meta
* First public release of the `mcradds` package.
* Submission to `CRAN`.

### New features
* Added `autoplot` method for Bland-Altman and regression plots.


# mcradds 0.2.0

### New features
* Added `tukey_outlier` and `dixon_outlier` to detect outliers ahead of establishing reference range.
* Added `robustRI` and `nonparRI` to compute robust and non-parametric reference range, and integrated into the main program `refInterval`.
* Wrapped `anovaVCA` and `VCAinference` from `VCA` package to analyze the variance components with ANOVA model.
* Added `aucTest` to do the AUC test for paired two-sample measurements in the designs of difference, non-inferiority and superiority.
* Added `RefInt` and `tpROC` classes and corresponding `show` method.
* Added `calcium`, `glucose`, `ldlroc` and `PDL1RP` data sets for example and testing use, and `nonparRanks` data set for internal function use.

### Enhancements
* Enhanced `diagTab` and `getAccuracy` so that they can support the reader precision ananlysis not only qualitative performance.

### Miscellaneous
* Added a series of helper function to format and concatenate to string.
* Uniform the capital and lower-case letters in roxygen documents.


# mcradds 0.1.0
* First release of the `mcradds` package, contains basic quantitative or qualitative performance methods and functions as shown below.

#### Sample Size
* Added `size_one_prop` and `size_ci_one_prop` for sample size of qualitative trials, `size_corr` and `size_ci_corr` for quantitative trails.

#### Classes and Datasets
* Added `SampleSize`, `MCTab` and `BAsummary` classes for `show` method.
* Added `platelet` and `qualData` data sets for example and testing use.

#### Analyzing Functions and Methods
* Added `diagTab` function to get 2x2 contingency table, and `getAccuracy` method to compute the qualitative diagnostic accuracy criteria.
* Added `blandAltman` function to calculate statistics for Bland-Altman, and `getOutlier` method to detect the potential outliers.
* Added `pearsonTest` and `spearmanTest`, more efficient functions to compute the confidence interval and hypothesis test.
* Added `mcreg` and `calcBias` methods of `mcr` package are wrapped for regression analysis.