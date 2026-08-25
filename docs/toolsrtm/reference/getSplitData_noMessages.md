# getSplitData for ML models with no messages

getSplitData for ML models with no messages

## Usage

``` r
getSplitData_noMessages(
  data = NULL,
  depVar = "Cab",
  inputs = NULL,
  data.trans = NULL,
  method.preProcess = NULL,
  depVar.trans = NULL,
  prop.split = NULL
)
```

## Arguments

- data:

  A dataframe with inputs and variable to predic

- depVar:

  name of the variable to predict

- inputs:

  the a vector with variables names for training and testing the model

- data.trans:

  a data.transformation method, options are: 'PCA','preProcess',

- method.preProcess:

  the data.transformation method for preProcess, data.transformation
  are: 'Normalize', 'YeoJohnson','BoxCox', Standarize',
  'Center','Scale', and 'PCA'

- depVar.trans:

  a boolean variable for applying data transformation in Y variable,
  options are: TRUE or FALSE.

- prop.split:

  a vector with proportion for spliting the dataset. prop.split
  =c(0.8,02) will be used as default.

## Value

the training and testing dataset with plot and scalers
