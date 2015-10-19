---
title: "Examples"
output: html_document
---

Below are some examples of how to use this package. We'll start be loading some datasets. The package comes with several pre-loaded datasets in scalar, interval, and monte carlo format. The hoverflies dataset (discussed by Spencer in his paper)

```{r}
summary(cars)
```

You can also embed plots, for example:

```{r, echo=FALSE}
plot(cars)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.


```{r, include=FALSE}
   # add this chunk to end of mycode.rmd
   file.rename(from="Examples.Rmd",to="README.md")
```
