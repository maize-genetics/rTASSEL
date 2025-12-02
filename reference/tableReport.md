# Return selected table report

Returns a `data.frame` object of association table reports

## Usage

``` r
tableReport(assocRes = missing(), reportName = missing())

# S4 method for class 'AssociationResults'
tableReport(assocRes = missing(), reportName = missing())

# S4 method for class 'AssociationResultsBLUE'
tableReport(assocRes = missing(), reportName = missing())

# S4 method for class 'AssociationResultsGLM'
tableReport(assocRes = missing(), reportName = missing())

# S4 method for class 'AssociationResultsMLM'
tableReport(assocRes = missing(), reportName = missing())

# S4 method for class 'AssociationResultsFast'
tableReport(assocRes = missing(), reportName = missing())

# S4 method for class 'AssociationResultsStepwise'
tableReport(assocRes = missing(), reportName = missing())

# S4 method for class 'PCAResults'
tableReport(assocRes = missing(), reportName = missing())
```

## Arguments

- assocRes:

  a
  [`AssociationResults`](https://rtassel.maizegenetics.net/reference/AssociationResults-class.md)
  object

- reportName:

  a specific table report to return
