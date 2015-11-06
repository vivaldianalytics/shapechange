---
title: "Function Descriptions and Examples"
output: html_document
---

Function: shape_change
----------------------

shape_change is the primary function in the package, and calculates any of 5 measures of shape change on 3 different types of data (scalar, interval, and monte carlo replicates). It takes 9 arguments:

* **data**: *Array*
<br/>
This argument should be a data frame, matrix, or 3D array. Each column contains abundances for each                         species, each row represents one time point, and (if a 3D array) the third dimension contains either two values (lower and upper bound for interval data) or > 2 values for Monte Carlo replicates.

* **time_val**: *String*
<br/>
This argument should be the column name of the dataset that represents time.

* **method**: *String*
<br/>
This tells the function which method to use when calculating shape change. Possible values:

    * Spencer
    * Foster and Tilman
    * Jassby and Goldman
    * Lewis
    * Bray Curtis Field
<br/>
* **surr**: *Boolean*
<br/>
If the data is scalar and the method chosen is *Spencer*, setting **sur** = TRUE will implement the surreal number calculations detailed by Spencer in the paper.

    
* **iter**: *Integer*
<br/>
If the data is in the form of Monte Carlo replicates, **iter** will determine how many of the replicates to run the calculations on. It is bounded below by 2 and above by the number of replicates contained in the data.

* **log_transform**: *Boolean*
<br/>
This variable tells the function whether or not the data is in log (base e) form.

* **mc_progress**: *Boolean*
<br/>
If the data is in the form of Monte Carlo replicates, this variable tells the function whether or not to display a progress bar as the simulations are run.

* **mc_int**: *Boolean*
<br/>
If the data is in the form of Monte Carlo replicates, this tells the function whether or not to return the results of calculations for all of replicates, or simply a Highest Posterior Density (HPD) confidence interval of the results (in which case the values returned would be intervals).

* **mc_prob**: *Double between 0 and 1*
<br/>
If **mc_int** == TRUE, this tells the function what the percent confidence of the HPD should be. The default value is 0.95, or 95% confidence.

Examples
--------

Below are some examples of how to use this function. We'll start be loading some datasets. The package comes with several pre-loaded datasets in scalar, interval, and monte carlo format. The hoverflies dataset (discussed by Spencer in his paper) contains abundance data for 14 species over a 30-year timeframe:


```r
data(hoverflies)
```

We can now calculate Spencer's shape change measure for the community using the shape_change function:


```r
library(shapechange)
SC_spencer = shape_change(data = hoverflies,"year",method = "Spencer")
```

On the other hand, let's suppose that we're dealing with interval data. In this case, we load the interval hoverflies dataset:


```r
data(hoverflies_INT)

SC_INT_spencer = shape_change(data = hoverflies_INT,"year",method = "Spencer")
```

Finally, we implement the function on the Monte Carlo dataset:


```r
data(hoverflies_MCMC)
SC_MCMC_spencer = shape_change(data = hoverflies_MCMC,"year",method = "Spencer")
```

```
##   |                                                                         |                                                                 |   0%
```


