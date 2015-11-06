


=======
# *shapechange* R Package

shapechange is an R package created to allow users to recreate the results obtained by Matthew M. Spencer in his paper titled 
*Size change, Shape change, and the growth space of a community*<sup>1</sup>. Spencer's paper 
proposed a new measure for the rate of shape change of a community. Using this package, users can calculate measures of change in size 
and shape of a community (respectively, proportional to among-species mean and standard deviation of mean proportional 
growth rates over some time interval). Size change is closely related to the Living Planet Index 
(see paper by Loh et al. 2005)<sup>2</sup>. In addition, users can calculate and plot other measures of change 
in communities, such as those proposed by Foster and Tilman<sup>3</sup>, Jassby and Goldman<sup>4</sup>, and Bray-Curtis<sup>5</sup>.





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
## 
  |                                                                       
  |                                                                 |   0%
```




# Citations

1: *Spencer, M. (2015). Size change, shape change, and the growth space of a community. Journal of Theoretical Biology 369:23-41. http://dx.doi.org/10.1016/j.jtbi.2015.01.002*

2: *Loh, J., Green, R. E., Ricketts, T., Lamoreux, J., Jenkins, M., Kapos, V., Randers, J., 2005. The Living Planet Index: using species population time series to track trends in biodiversity. Philosophical Transactions of the Royal Society B 360, 289–295.*

3: *Foster, B. L., Tilman, D., 2000. Dynamic and static views of succession: Testing the descriptive power of the chronosequence approach. Plant Ecology 146, 1–10.*

4: *Jassby, A. D., Goldman, C. R., 1974. A quantitative measure of succession rate and its application to the phyto- plankton of lakes. American Naturalist 108, 688–693.*

5: *Field, J. G., Clarke, K. R., Warwick, R. M., 1982. A practical strategy for analysing multispecies distribution patterns. Marine Ecology Progress Series 8, 37–52.*
>>>>>>> 569e86f62149defbf1636f407bec820ee52dfc7c
