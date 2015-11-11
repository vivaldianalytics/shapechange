You can install the library straight from gihub using the following lines of code. Installation requires the package devtools.

*require(devtools)* <br/> *\# If devtools is not installed, enter: install.packages("devtools")* <br/> *install\_github("vivaldianalytics/shapechange")*

Function Descriptions and Examples
----------------------------------

*shape\_change*
---------------

shape\_change is the primary function in the package, and calculates any of 5 measures of shape change on 3 different types of data (scalar, interval, and monte carlo replicates). It takes 9 arguments:

-   **data**: *Array* <br/> This argument should be a data frame, matrix, or 3D array. Each column contains abundances for each species, each row represents one time point, and (if a 3D array) the third dimension contains either two values (lower and upper bound for interval data) or \> 2 values for Monte Carlo replicates.

-   **time\_val**: *String* <br/> This argument should be the column name of the dataset that represents time.

-   **method**: *String* <br/> This tells the function which method to use when calculating shape change. Possible values:

    -   Spencer
    -   Foster and Tilman
    -   Jassby and Goldman
    -   Lewis
    -   Bray Curtis Field <br/>
-   **surr**: *Boolean* <br/> If the data is scalar and the method chosen is *Spencer*, setting **sur** = TRUE will implement the surreal number calculations detailed by Spencer in the paper.

-   **iter**: *Integer* <br/> If the data is in the form of Monte Carlo replicates, **iter** will determine how many of the replicates to run the calculations on. It is bounded below by 2 and above by the number of replicates contained in the data.

-   **log\_transform**: *Boolean* <br/> This variable tells the function whether or not the data is in log (base e) form.

-   **mc\_progress**: *Boolean* <br/> If the data is in the form of Monte Carlo replicates, this variable tells the function whether or not to display a progress bar as the simulations are run.

-   **mc\_int**: *Boolean* <br/> If the data is in the form of Monte Carlo replicates, this tells the function whether or not to return the results of calculations for all of replicates, or simply a Highest Posterior Density (HPD) confidence interval of the results (in which case the values returned would be intervals).

-   **mc\_prob**: *Double between 0 and 1* <br/> If **mc\_int** == TRUE, this tells the function what the percent confidence of the HPD should be. The default value is 0.95, or 95% confidence.

Examples
--------

Below are some examples of how to use this function. We'll start be loading some datasets. The package comes with several pre-loaded datasets in scalar, interval, and monte carlo format. The hoverflies dataset (discussed by Spencer in his paper) contains abundance data for 14 species over a 30-year timeframe:

``` {.r}
data(hoverflies)
```

We can now calculate Spencer's shape change measure for the community using the shape\_change function:

``` {.r}
library(shapechange)
SC_spencer = shape_change(data = hoverflies,"year",method = "Spencer")
```

On the other hand, let's suppose that we're dealing with interval data. In this case, we load the interval hoverflies dataset:

``` {.r}
data(hoverflies_INT)

SC_INT_spencer = shape_change(data = hoverflies_INT,"year",method = "Spencer")
```

Finally, we implement the function on the Monte Carlo dataset:

``` {.r}
data(hoverflies_MCMC)
SC_MCMC_spencer = shape_change(data = hoverflies_MCMC,"year",method = "Spencer")
```

    ## 
      |                                                                       
      |                                                                 |   0%

*as.interval.frame*
-------------------

as.interval.frame converts scalar datasets into interval-valued datasets using one of two methods: Square Root Method (send each x to x +/- p\*sqrt(x)), and the Poisson method (calculating a poisson confidence interval about x). This function takes 4 arguments:

-   **data**: *Data Frame or Matrix* <br/> This argument should be a data frame or a matrix. Each column contains abundances for each species, each row represents one time point. An additional column denoting time is also allowed, however the function should be informed as to which column it is:

-   **exclude**: *String* <br/> This argument tells the function which column to exclude from the transformation process. Usually this column should hold values denoting time.

-   **method**: *String* <br/> Tells the function which method to use (Square Root or Poisson) to calculate the confidence intervals. Possible values:

-   'sqrt'
-   'poisson'

-   **param**: *double* <br/> Both of the previous methods described take one parameter. The Square Root method, which sends each point x to x +/- p\*sqrt(x), typically should have a parameter p of around 2. The Poisson method calculates a confidence interval based on a percentage confidence. If the user would like the function to return 95% confidence intervals, they should set the param value to 0.95.

Examples
--------

Here are a couple of demonstrations of how to use this function.

Using the square root method:

``` {.r}
hoverflies_sqrt = as.interval.frame(hoverflies,exclude = "year",method = "sqrt",param = 2)
```

Using the Poisson method:

``` {.r}
hoverflies_pois = as.interval.frame(hoverflies,exclude = "year",method = "poisson",param = 0.95)
```
