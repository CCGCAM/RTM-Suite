# getMLmodel with a prefixed configuration

getMLmodel with a prefixed configuration

## Usage

``` r
getMLmodel(
  dataset = NULL,
  depVar = "Cab",
  model = "CNN",
  optimizer = "adam",
  batch.size = 125,
  n.epochs = 10,
  save.model = T,
  path.model = NULL,
  prop.split = c(0.8, 0.2),
  data.trans = "preProcess",
  method.preProcess = "Normalize",
  depVar.trans = FALSE
)
```

## Arguments

- dataset:

  a dataframe

- depVar:

  name of the variable to predict

- model:

  a ML model. options are: 'CNN','Hidden-layers',

- optimizer:

  the optimizer for the model. options are: 'adam','adadelta','adagrad',
  'adamax', 'nadam', 'msprop', 'sgd'

- batch.size:

  batch size used for each epoch. By default is 125

- n.epochs:

  number of epoch. By default is 100

- save.model:

  a boolean variable for saving ML model, options are: TRUE or FALSE. if
  TRUE, please use path.model to give a folder for the model

- path.model:

  a path for saving the models. By default path.model ='Models'

- prop.split:

  a vector with proportion for spliting the dataset. prop.split
  =c(0.8,02) will be used as default.

- data.trans:

  a data.transformation method, options are: 'PCA','preProcess',

- method.preProcess:

  the data.transformation method for preProcess, data.transformation
  are: 'Normalize', 'YeoJohnson','BoxCox', Standarize',
  'Center','Scale', and 'PCA'

- depVar.trans:

  a boolean variable for applying data transformation in Y variable,
  options are: TRUE or FALSE.

## Value

a list with models and plots
