

rm(list= ls())
# 1. load the main libraries  -----

library(ToolsRTM)
required_packages <- c("shiny", "shinythemes", 'shinydashboard',"ggplot2", "dplyr", "tidyr",
                       'shinybusy','parallel',"doParallel",'foreach',
                       "tensorflow", "keras", "caret", "e1071", "viridis")


# Check for missing packages and install them if necessary
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}
# Load the libraries
lapply(required_packages, library, character.only = TRUE)


# 2. Define UI -----

ui <- navbarPage("Online reflectance simulator",theme = shinytheme("flatly"),
                 
                 # 5) Tab Inversion module  -----------------------------------------
                 
                 tabPanel(title = "Inversion of Plant Traits via Hybrid ML models",
                      
                          # Add a sidebar layout
                          sidebarLayout(
                            # Add a sidebar panel
                            sidebarPanel(
                              # Add information on cranes and prompt user to use slider
                              h4("Select the parameters:"),
                              
                              # LUT database selection
                              selectInput("lut_db", "Select LUT Database:",
                                          choices = list(
                                            "Upload New Dataset" = "upload",
                                            "PROSAILH (pre-defined LUT)" = "PROSAIL",
                                            "INFORM (pre-defined LUT)" = "INFORM"
                                          )),
                              # Conditional panel to show file input when "Upload New Dataset" is selected
                              conditionalPanel(
                                condition = "input.lut_db == 'upload'",
                                
                                fileInput("uploaded_file", "Upload Your Dataset (CSV):",
                                          accept = c(".csv"), multiple = FALSE)
                                
                              ),
                              # Checkbox for plotting indices
                              checkboxInput("plotting_indices", "Calculate spectral indices", value = FALSE),
                              
                              
                              tags$hr(style = "border-color: grey;"),
                              
                              #h4("Invget.saersion Method:"),
                              # Dependent variable selection
                              selectInput("depVar", "Select a plant trait to estimate:",
                                          choices = c('Chorophyll content (Cab)' ='Cab',
                                                      'Carotenoid content (Car)' ='Car', 'Anthocyanin  content (Anth)' ='Anth',
                                                      'Leaf area index (LAI)'= 'LAI', 'leaf water content (EWT)' = 'EWT',
                                                      'leaf protein content (Prot)' = 'Prot',
                                                      'Carbon-based constituent (CBC)' = 'CBC',
                                                      'Leaf brown pigments (Cbrown)'= 'Cbrown')),
                              
                              tags$hr(style = "border-color: grey;"),
                              
                              selectInput("inputs_toML", "Select input features for training the ML model:",
                                          choices = c( 'Red-edge Vegetation Indices' = 'rededge_indices',
                                                       'Main Vegetation Indices' = 'indices',
                                                       'Visble-NIR Vegetation Indices' = 'visible_indices',
                                                       
                                                       'SWIR Vegetation Indices' = 'swir_indices',
                                                       'Reflectance Bands' = 'reflectance',
                                                       'Reflectance VNIR Bands' = 'vnir_reflectance',
                                                       'Reflectance SWIR Bands' = 'swir_reflectance',
                                                       'Combined Vegetation Indices and Reflectance' = 'indices_reflectance',
                                                       'Combined Vegetation Indices and Reflectance in VNIR' = 'indices_reflectance_vnir',
                                                       'Combined Vegetation Indices and Reflectance in SWIR' = 'indices_reflectance_swir')),
                              
                              
                              tags$hr(style = "border-color: grey;"),
                              
                              
                              
                              # Inversion model selection
                              selectInput("models_", "Select Model:",
                                          choices = list(
                                            "Support Vector Machine" = "SVM",
                                            "deep Neural Network" = "Hidden_layers",
                                            "Random Forest" = "RF",
                                            #        "Convolution Neural Network" = "CNN",
                                            #    "Neural Network (caret)" = "NN",
                                            "Partial least square regression" = "PLSR",
                                            "eXtreme Gradient Boosting (XGBoost)" = "xGB",
                                            "Gradient Boosting" = "GB"
                                          )),
                              
                              
                              
                              
                              tags$hr(style = "border-color: grey;"),
                              
                              # Conditional panel for neural networks (CNN or NNe)
                              conditionalPanel(
                                condition = "input.models_ == 'CNN' || input.models_ == 'Hidden_layers'",
                                
                                h4("Neural Network Configuration"),
                                
                                # Preprocessing method selection
                                selectInput("method.preProcess", "Select Preprocessing Method:",
                                            choices = c('Normalize', 'Standarize', 'Center',
                                                        'YeoJohnson', 'BoxCox')),
                                
                                # Optimizer selection
                                selectInput("optimizer", "Select Optimizer:",
                                            choices = c('adam', 'adadelta', 'adagrad', 'adamax',
                                                        'nadam', 'msprop', 'sgd')),
                                
                                # Neural network parameters
                                numericInput("n_layers", "Number of Layers:", 3, min = 1,max=5),
                                numericInput("n_neurons", "Number of Neurons:", 128, min = 1,max=1024),
                                numericInput("batch_size", "Batch Size:", 32, min = 1,max=2048),
                                numericInput("n_epochs", "Number of Epochs:", 10, min = 1, max=1000),
                                numericInput("p_samplesML_keras", "Percentage of samples:", value=1, min = 1,max=90),
                              ), # Conditional panel for neural networks (CNN or NNe)
                              conditionalPanel(
                                condition = "input.models_ == 'SVM' || input.models_ == 'RF'  || input.models_ == 'NN' || input.models_ == 'PLSR'
                                || input.models_ == 'xGB' || input.models_ == 'GB'",
                                
                                h4("Parameters: "),
                                
                                # Preprocessing method selection
                                numericInput("p_samplesML", "Percentage of samples:", value=1, min = 1,max=90),
                                p(em('Adjust this value to control the proportion of the samples used in the training processing.
                                     Increasing the percentage may enhance the accuracy of the results but will require more processing time.'))
                                
                              ),
                              
                              actionButton("train_model", "Train Model")
                              
                              
                              
                            ), # Close sidebarPanel
                            mainPanel( # open main panel
                              mainPanel(
                                tabsetPanel(
                                  id = "tabs",
                                  
                                  tabPanel("LUTs",
                                           h4(""),
                                           # Conditional panel for showing the outputs based on the selected LUT database
                                           conditionalPanel(
                                             condition = "input.lut_db != ''",  # Trigger panel when any LUT database is selected
                                             h4('Statistical Reflectance Profiles for the ML Inversion'),
                                             # Display a custom error message if there's no uploaded file or processing issues
                                             uiOutput("error_message"), # Placeholder for error message
                                             
                                             # Output for the plot
                                             plotOutput("plot_checking"),
                                             br(),
                                             # uiOutput("lut_selection"),  # Links based on `notebook_p`
                                             
                                             # Output for the DataTable
                                             DT::dataTableOutput("lut_header"),
                                           ),
                                           
                                  ), ## end panel
                                  
                                  tabPanel("Indices",
                                           h4("Scatterplot Analysis"),
                                           # Informative note for users
                                           p("Note: Estimation of indices is required. Please ensure to check the 'Calculate Indices' box to enable this feature."),
                                           br(),
                                           # UI for selecting variables for the scatter plot
                                           uiOutput("select_axes_ui"),  # Dynamically generated UI for axis selection
                                           
                                           # Output for the scatter plot
                                           plotOutput("scatter_plot"),  # Render the scatter plot
                                           br(),
                                           # Informative note for users
                                           uiOutput("instructions_ui"),  # Dynamically generated UI for instructions
                                           
                                           
                                           
                                           
                                  ), ### close tab panel
                                  tabPanel("Predictions",
                                           h3(""),
                                           br(),
                                           p('This section displays a scatterplot of predicted
                                             inputs (e.g., Cab) based on the selected ML model.
                                             The plot visualizes model predictions using the selected inputs from the
                                             testing dataset, along with key statistics for performance evaluation.'),
                                           
                                           
                                           textOutput("status"),
                                           verbatimTextOutput("training_progress"),  # Real-time training progress
                                           DT::DTOutput("training_history") , # Final training history
                                           # DT::dataTableOutput("training_history"),  # New output for training history
                                           plotOutput("predictionPlot"),
                                           br(),
                                           DT::dataTableOutput("statsTable_ML"),
                                           br(),
                                           uiOutput("plotReady_tableReady"),
                                           
                                           br(),
                                           tags$hr(style = "border-color: grey;"),
                                           
                                           
                                  ), ## end panel
                                  
                                  tabPanel("Info",
                                           
                                           h4("How use the inversion module"),
                                           # Custom CSS for instruction box
                                           tags$style(HTML("
                                           .instruction-box {
                                           border: 1px solid #007BFF;  /* Blue border */
                                           background-color: #F0F8FF;  /* Light blue background */
                                           padding: 15px;
                                           margin-bottom: 20px;
                                           border-radius: 5px;
                                           font-size: 16px;
                                           color: #333;
                                           }
                                                     ")),
                                           
                                           # Instruction Box
                                           p(strong("Step 1:")),
                                           
                                           p(style = "text-align: justify;",
                                             "Select a pre-trained Look-up table (LUT) for uploading your dataset in the required format accroding the LUT module."),
                                           
                                           p(em("Note₁: Review the spectral profiles of your reflectance bands in the 'LUTs' tab.")),
                                           
                                           # Instruction Box
                                           p(strong("Step 2:")),
                                           # Instruction Box
                                           
                                           p(style = "text-align: justify;",
                                             "Estimate a set of key spectral indices related to biochemical, biophysical, and structural
                                                      traits for Sentinel-2 resolution data"),
                                           
                                           p(em("Note₂: Check the relationship between indices and model inputs in the 'Indices' tab.")),
                                           p(em('Note₃: These indices will be used to train a model for retrieving key plant traits.')),
                                           
                                           # Instruction Box
                                           p(strong("Step-3:")),
                                           
                                           p(style = "text-align: justify;",
                                             "Select a deep machine learning (ML) or classical ML model and adjust the parameters."),
                                           
                                           p(em("Note₄: Once the deep learning model parameters are
                                           defined, proceed to train the model by clicking the 'Train the model' button.")),
                                           p(em("Note₅: The deep ML algorithms are the 1D-Convolutional Neural Networks (CNN) and a clasical ML model with hidden layers.")),
                                           
                                           p(strong("Step-4:")),
                                           
                                           p(style = "text-align: justify;",
                                             "In the 'Predictions' tab, review the results and download the scaler, statistical scores, and the trained model,
                                                     "),
                                           
                                           p(em("This module utilizes the LUT generated by the PROSAIL-PRO model (PROSPECT-PRO & fourSAILH) as described in ", strong("Camino et al. (2024, in preparation)."))),
                                           br(),
                                           # Instruction Box
                                           tags$div(class = "instruction-box",
                                                    p(style = "text-align: justify; font-size: 12px;",
                                                      "This module is designed to work with Sentinel-2A and Sentinel-2B imagery using a convolved LUT. Herem thw plant traits were retrieved from reflectance data
                                                     ad main spectral indices followed the method described in Camino et al., (2024; in prep.)
                                                     entitled 'Integrating physiological plant traits with Sentinel-2 imagery for monitoring gross primary production and detecting forest disturbances.'"),
                                                    p(style = "text-align: justify; font-size: 12px;",
                                                      "Future updates will integrate the ", strong("ToolsRTMs")," package, enabling compatibility with additional platforms and expanding remote sensing applications. (e.g., MODIS, Landsat)"),
                                                    p(style = "text-align: justify; font-size: 12px;",
                                                      "For support with other satellite (e.g., Landsat) or aerial platforms (e.g., hyperspectral imagery), please consult the ", strong("ToolsRTMs")," package documentation on CRAN.")
                                           ),
                                           br(),
                                           p(em("The module uses the ", code("getMLmodel.withRetrain")," function, integrated into ", strong("ToolsRTMs")," package.")),
                                           p(em("The module uses the ", code("get.inversion")," function, integrated into ", strong("ToolsRTMs")," package.")),
                                           
                                           
                                  ), ## end panel
                                  
                                  tabPanel("Functions",
                                           # Introduction to the module
                                           h4("Main funtions integrated in the ToolsRTM Package"),
                                           
                                           p(style = "text-align: justify;", "This module utilizes key functions integrated into the ", strong("ToolsRTM"), " package to process and analyze vegetation data derived from Sentinel-2 satellite imagery. These functions are critical for estimating spectral indices, training machine learning models, and predicting plant traits such as gross primary production (GPP) and other physiological traits."),
                                           br(),
                                           # Section for getIndicesSE2
                                           h4("1. getIndicesSE2"),
                                           p(style = "text-align: justify;", "The ", code("getIndicesSE2"), " function is designed to compute various spectral indices using Sentinel-2 data. It processes the hyperspectral bands to derive essential vegetation indices like ", strong("NDVI"), ", ", strong("EVI"), ", and others. These indices are crucial for estimating parameters such as GPP and monitoring the health and traits of vegetation."),
                                           
                                           # Arguments
                                           p(strong("KeyParameters:")),
                                           tags$ul(
                                             tags$li(code("df"), " - A dataframe with reflectance where each row corresponds to a spectrum."),
                                             tags$li(code("sensor"), " - Sensor options: 'Sentinel-2a' or 'Sentinel-2b'."),
                                             tags$li(code("df.data"), " - Dataset with IDs corresponding to each spectrum; can be NULL."),
                                             tags$li(code("fast.process"), " - Set to TRUE if bands are ordered for Sentinel-2; otherwise, FALSE or NULL.")
                                           ),
                                           # Usage
                                           p(strong("Usage:")),
                                           code("getIndicesSE2(df, sensor = 'Sentinel-2a', df.data = NULL, "),
                                           br(),
                                           code("fast.process = NULL)"),
                                           br(),
                                           br(),
                                           # Section for getMLmodel.withRetrain
                                           h4("2. getMLmodel.withRetrain"),
                                           p(style = "text-align: justify;", "The ", code("getMLmodel.withRetrain"), " function is a deep learning-based model retraining process. It helps fine-tune machine learning models by retraining them on the input features and vegetation indices derived from Sentinel-2 data. The model is trained to predict variables like GPP or specific plant traits. This function depends on the following packages: ", code("MLMetric"), ", ", code("keras"), ", and ", code("tensorflow"), ". Note that the non-retraining version of this function is ", code("getMLmodel"), " from the ToolsRTM package."),
                                           
                                           # Key parameters explanation
                                           p(strong("Key Parameters:")),
                                           p(style = "text-align: justify;", code("dataset"), " - The dataset that contains the input features (bands) and the dependent variable."),
                                           p(style = "text-align: justify;", code("depVar"), " - The dependent variable (e.g., GPP, LAI, chlorophyll content) being predicted."),
                                           p(style = "text-align: justify;", code("model"), " - The type of machine learning model to be used. Options include 'CNN' for Convolutional Neural Networks or 'Hidden-layers' for models with specified hidden layers."),
                                           p(style = "text-align: justify;", code("optimizer"), " - The optimization algorithm used for training the model. Options include 'adam', 'adadelta', 'adagrad', 'adamax', 'nadam', 'msprop', and 'sgd'."),
                                           p(style = "text-align: justify;", code("n.times"), " - The number of times the model retrains. By default, this is set to 1, but it can be adjusted for multiple training iterations."),
                                           p(style = "text-align: justify;", code("n.neurons"), " - The number of neurons to use in the hidden layers of the model. Default is 128 neurons, applicable primarily for models that utilize hidden layers."),
                                           p(style = "text-align: justify;", code("n.layers"), " - The number of hidden layers in the model. Default is 4, only relevant for the 'Hidden-layers' model."),
                                           p(style = "text-align: justify;", code("batch.size"), " - The number of samples processed before the model is updated. Default is 125, but it can be adjusted based on the dataset size."),
                                           p(style = "text-align: justify;", code("n.epochs"), " - The number of complete passes through the training dataset. Default is 100 epochs, which can be modified as needed."),
                                           p(style = "text-align: justify;", code("save.model"), " - A boolean value that specifies whether to save the trained model. Options are TRUE or FALSE; if TRUE, specify a path in ", code("path.model"), " for saving."),
                                           p(style = "text-align: justify;", code("path.model"), " - The path for saving the trained model. If not specified, defaults to 'Models'."),
                                           p(style = "text-align: justify;", code("prop.split"), " - A vector indicating the proportion of the dataset to be used for training and validation. The default is ", code("c(0.8, 0.2)"), " which represents 80% for training and 20% for validation."),
                                           p(style = "text-align: justify;", code("data.trans"), " - Specifies the method of data transformation to be applied to the dataset before training. Options include 'PCA' and 'preProcess'."),
                                           p(style = "text-align: justify;", code("method.preProcess"), " - The specific data transformation method applied during preprocessing. Options include 'Normalize', 'YeoJohnson', 'BoxCox', 'Standardize', 'Center', 'Scale', and 'PCA'."),
                                           p(style = "text-align: justify;", code("depVar.trans"), " - A boolean value indicating whether to apply data transformation to the dependent variable. Options are TRUE or FALSE."),
                                           br(),
                                           p(strong("Return Value:")),
                                           p(style = "text-align: justify;", "The ", code("getMLmodel.withRetrain"), " function returns a list with the following elements:"),
                                           tags$ul(
                                             tags$li(code("model"), " - The trained Keras model object."),
                                             tags$li(code("history"), " - The training history of the model (loss and accuracy metrics over epochs)."),
                                             tags$li(code("stats"), " - Performance statistics to export (e.g., RMSE, R²)."),
                                             tags$li(code("Scalar.train"), " - The scaling parameters applied to the training data."),
                                             tags$li(code("Scalar.Ytrain"), " - The scaling parameters applied to the dependent variable for the training data."),
                                             tags$li(code("plot.val"), " - A plot visualizing validation results."),
                                             tags$li(code("plot.cor"), " - A plot displaying correlation between predicted and actual values.")
                                           ),
                                           p(style = "text-align: justify;", "An example of the model's performance visualization is shown below:"),
                                           img(src = "images/model_cab.png", alt = "Model RandomForest Performance", style = "width: 100%; height: auto;"),
                                           br(),
                                           br(),
                                           p(strong("Usage:")),
                                           # Example usage of the getMLmodel.withRetrain function
                                           p(style = "text-align: justify;", "Example of using the ", code("getMLmodel.withRetrain"), " function:"),
                                           code(style = "text-align: justify;","trained_model <- getMLmodel.withRetrain("),
                                           code(style = "text-align: justify;","dataset = my_data_frame, depVar = 'Cab', model = 'CNN', optimizer = 'adam',"),
                                           code(style = "text-align: justify;","n.times = 1, n.neurons = 128, n.layers = 3, batch.size = 32,"),
                                           code(style = "text-align: justify;","n.epochs = 100, save.model = TRUE, path.model = 'Models/',"),
                                           code(style = "text-align: justify;","prop.split = c(0.8, 0.2), data.trans = 'preProcess',"),
                                           code("method.preProcess = 'Normalize', depVar.trans = FALSE)"),
                                           br(),
                                           # Section for get.inversion
                                           br(),
                                           h4("3. get.inversion"),
                                           p(style = "text-align: justify;", "The ", code("get.inversion"), " function performs inversion of plant traits using various machine learning models. The function relies on several packages, including ", code("caret"), ", ", code("randomForest"), ", ", code("e1071"), ", ", code("doParallel"), ", ", code("caretEnsemble"), ", and ", code("stats"), ". It allows you to run different machine learning algorithms for predicting dependent variables like plant traits from independent variables (inputs). Below are the key configurations for this function:"),
                                           
                                           # Key parameters explanation
                                           p(strong("Key Parameters:")),
                                           p(style = "text-align: justify;", code("data"), " - The dataset containing the independent and dependent variables."),
                                           p(style = "text-align: justify;", code("depVar"), " - The dependent variable that you want to predict (e.g., GPP, LAI, chlorophyll content)."),
                                           p(style = "text-align: justify;", code("inputs"), " - A vector of independent variables (e.g., vegetation indices, reflectance bands) used for training the model."),
                                           p(style = "text-align: justify;", code("algorithm"), " - The type of machine learning model. Available options include:"),
                                           tags$ul(
                                             tags$li(code("PLSR"), " - Partial Least Squares Regression"),
                                             tags$li(code("SVM"), " - Support Vector Machine"),
                                             tags$li(code("RF"), " - Random Forest"),
                                             tags$li(code("NN"), " - Neural Network"),
                                             tags$li(code("GB"), " - Gradient Boosting"),
                                             tags$li(code("xGB"), " - Extreme Gradient Boosting (XGBoost) with linear base learners"),
                                             tags$li(code("Bayesian"), " - Bayesian Additive Regression Trees"),
                                             tags$li(code("AdaBag"), " - Bagged AdaBoost"),
                                             tags$li(code("qLASSO"), " - Quantile Regression with LASSO penalty"),
                                             tags$li(code("RVM"), " - Relevance Vector Machines (RVM) with linear kernel"),
                                             tags$li(code("BRNN"), " - Bayesian Regularized Neural Networks"),
                                             tags$li(code("Ensemble"), " - Stacking Ensemble models")
                                           ),
                                           p(style = "text-align: justify;", code("method.resampling"), " - The resampling method for controlling the training of the ML model. Options include 'boot' (Bootstrapping), 'cv' (Cross-Validation), 'LOOCV' (Leave-One-Out Cross-Validation), among others."),
                                           p(style = "text-align: justify;", code("n.cores"), " - The number of CPU cores to use for parallel processing. If not specified, defaults to single-core."),
                                           p(style = "text-align: justify;", code("seed"), " - The seed for random number generation, used to ensure reproducibility of results. Default is 123."),
                                           p(style = "text-align: justify;", code("n.samples"), " - The number of samples used for tuning and model training. Default is 500."),
                                           p(style = "text-align: justify;", code("save.model"), " - A logical value indicating whether to save the trained model. If TRUE, specify the path using ", code("save.path"), ". Default is FALSE."),
                                           p(style = "text-align: justify;", code("save.path"), " - The file path where the trained model will be saved, if ", code("save.model"), " is set to TRUE."),
                                           
                                           br(),
                                           p(strong("Return Value:")),
                                           p(style = "text-align: justify;", "The ", code("get.inversion"), " function returns a list with the following elements:"),
                                           tags$ul(
                                             tags$li(code("model.label"), " - The name of the algorithm used (e.g., 'SVM', 'RF')."),
                                             tags$li(code("model"), " - The trained model object."),
                                             tags$li(code("predictions"), " - A list containing the model's predictions for both the training and testing datasets:"),
                                             tags$ul(
                                               tags$li(code("train"), " - Predictions for the training dataset."),
                                               tags$li(code("test"), " - Predictions for the testing dataset.")
                                             ),
                                             tags$li(code("statistics"), " - Model performance statistics (e.g., RMSE, R², etc.)."),
                                             tags$li(code("plot"), " - A plot showing model performance or predictions."),
                                             tags$li(code("importance"), " - Variable importance scores (if applicable, depending on the model type).")
                                           ),
                                           p(style = "text-align: justify;", "An example of the model's performance visualization is shown below:"),
                                           img(src = "images/model_xgb.png", alt = "Model XGBoost Performance", style = "width: 100%; height: auto;"),
                                           
                                           br(),
                                           p(strong("Usage:")),
                                           # Example usage of the get.inversion function
                                           p(style = "text-align: justify;", "Example of using the ", code("get.inversion"), " function:"),
                                           code(style = "text-align: justify;", "inversion_results <- get.inversion("),
                                           code(style = "text-align: justify;", "data = my_data, depVar = 'Cab', inputs = c('NDVI', 'TCARI'), algorithm = 'SVM',"),
                                           code(style = "text-align: justify;", "method.resampling = 'cv', seed = 123, n.samples = 500, save.model = TRUE, save.path = 'Models/')"),
                                           br(),
                                           br(),
                                           # Datasets used for Sentinel-2 applications
                                           h4("Datasets used for Sentinel-2 Applications"),
                                           p("The datasets utilized in these functions include both pre-established and specifically developed datasets for Sentinel-2 applications:"),
                                           
                                           # PROSAIL Dataset
                                           p(style = "text-align: justify;", strong("PROSAIL Dataset:"), " The PROSAIL dataset serves as the Look-Up Table (LUT) in the study ", em("Camino et al. 2024 (in prep.)"), " for integrating physiological plant traits with Sentinel-2 imagery. It is widely used for estimating variables like LAI, chlorophyll content, and other key plant traits."),
                                           
                                           # INFORM Dataset
                                           p(style = "text-align: justify;", strong("INFORM Dataset:"), " The INFORM dataset is developed for a study on bark beetle outbreaks in forest canopies. This dataset is tailored for analyzing the impact of bark beetles on forest health using Sentinel-2 imagery, and it forms the basis of a study currently in preparation (", em("Camino et al. 2025 (in prep.)"), ")."),
                                           # Section for SPART
                                           h4("Citation"),
                                           
                                           p(style = "text-align: justify;",HTML('If you use these functions with <b>ToolsRTM</b> package, please cite the following references:')),
                                           
                                           p(style = "text-align: justify;",HTML('Camino et al., (2024). RT-Simulator: An Online Platform to Simulate Canopy Reflectance from Biochemical and Structural Plant Properties Using Radiative Transfer Models,
                                          <i>IGARSS 2024 - 2024 IEEE International Geoscience and Remote Sensing Symposium</i>, Athens, Greece, 2024, pp. 2811-2814,
                                          <a href="https://ieeexplore.ieee.org/document/10642442" target="_blank">doi: 10.1109/IGARSS53475.2024.10642442</a>.')),
                                           
                                           p(style = "text-align: justify;",HTML('Arano et al., (2024). Enhancing Chlorophyll Content Estimation With Sentinel-2 Imagery: A Fusion of Deep Learning and Biophysical Models,
                                          <i>IGARSS 2024 - 2024 IEEE International Geoscience and Remote Sensing Symposium</i>, Athens, Greece, 2024, pp. 4486-4489,
                                          <a href="https://ieeexplore.ieee.org/document/10641613" target="_blank">doi: 10.1109/IGARSS53475.2024.10641613</a>.')),
                                           
                                           p(style = "text-align: justify;",'Camino et al., (in prep). Integrating physiological plant traits with
                                             Sentinel-2 imagery for monitoring gross primary production and detecting forest disturbances. '),
                                           
                                           br()
                                  ),
                                  tabPanel("References",
                                           h3(""),
                                           
                                           p(style = "text-align: justify;",HTML('If you use <b>ToolsRTM</b> or <b>SCOPEinR</b> packages, please cite the following references:')),
                                           
                                           p(style = "text-align: justify;",HTML('Camino et al., (2024). RT-Simulator: An Online Platform to Simulate Canopy Reflectance from Biochemical and Structural Plant Properties Using Radiative Transfer Models,
                                   <i>IGARSS 2024 - 2024 IEEE International Geoscience and Remote Sensing Symposium</i>, Athens, Greece, 2024, pp. 2811-2814,
                                          <a href="https://ieeexplore.ieee.org/document/10642442" target="_blank">doi: 10.1109/IGARSS53475.2024.10642442</a>.')),
                                           
                                           p(style = "text-align: justify;",HTML('Arano et al., (2024). Enhancing Chlorophyll Content Estimation With Sentinel-2 Imagery: A Fusion of Deep Learning and Biophysical Models,
                                   <i>IGARSS 2024 - 2024 IEEE International Geoscience and Remote Sensing Symposium</i>, Athens, Greece, 2024, pp. 4486-4489,
                                   <a href="https://ieeexplore.ieee.org/document/10641613" target="_blank">doi: 10.1109/IGARSS53475.2024.10641613</a>.')),
                                           
                                           p(style = "text-align: justify;",'Camino et al., (in prep). Integrating physiological plant traits with
                                             Sentinel-2 imagery for monitoring gross primary production and detecting forest disturbances. '),
                                           br(),
                                  ) ### close tab panel
                                )  ## close tabsetPanel
                                
                                
                              ) # Close mainPanel
                            ) # Close mainPanel
                          ) # Close sideBarLayout
                 ), # Close the breeding distribution tab panel
                 # Add a References tab panel, linking each Reference to its source.
                 
                 
                 
)
# 3. Define server logic -----


server <- function(input, output,session) {
  
  #Sys.setenv(TF_CPP_MIN_LOG_LEVEL = '2')  # Set log level to reduce TensorFlow output
  #Sys.setenv("CUDA_VISIBLE_DEVICES" = "-1")  # Disable GPU
  library(tensorflow)
  tf$config$optimizer$set_jit(FALSE)  # Disable XLA JIT compiler
  # 5) MACHINE LEARNIG module   ---------------------------------------
  
  
  # Reactive value to keep track of whether the dataset has been uploaded
  dataset_uploaded <- reactiveVal(FALSE)
  
  # Create a reactive value to store the uploaded dataset
  user_lut_data <- reactiveVal()
  # Create a reactive value to store the uploaded dataset
  user_lut_data_indices <- reactiveVal()
  
  ## Get datasets   ---------------------------------------
  
  
  loadDataset <- function(dataset_type) {
    if (dataset_type == "upload" && !is.null(input$uploaded_file)) {
      # Load uploaded dataset
      user_lut <- tryCatch({
        read.csv(input$uploaded_file$datapath)
      }, error = function(e) {
        showNotification("Error loading uploaded file: Please check the file format.", type = "error")
        return(NULL)
      })
      if (!is.null(user_lut)) {
        showNotification("File uploaded successfully.", type = "message")
      }
      return(user_lut)
      
    } else if (dataset_type == "PROSAIL") {
      # Load predefined PROSAIL dataset
      user_lut <- readRDS("www/datasets/LUT_prosailSE2.rds")
      showNotification("PROSAIL dataset loaded successfully.", type = "message")
      return(user_lut)
      
    } else if (dataset_type == "INFORM") {
      # Load predefined INFORM dataset
      user_lut <- readRDS("www/datasets/LUT_informSE2.rds")
      showNotification("INFORM dataset loaded successfully.", type = "message")
      return(user_lut)
      
    } else {
      # No dataset selected or uploaded
      showNotification("No valid dataset selected or uploaded.", type = "error")
      return(NULL)
    }
  }
  
  
  ## Get Processed datasets   ---------------------------------------
  
  # Reactive function to load the selected or uploaded dataset
  processed_data <- reactive({
    req(input$lut_db)  # Ensure `lut_db` input is available
    dataset <- loadDataset(input$lut_db)
    user_lut_data(dataset)
    return(dataset)
  })
  # UI Output for custom error message
  output$error_message <- renderUI({
    if (is.null(processed_data())) {
      div(style = "color: red; font-weight: bold;",
          "Error: No data available. Please upload a valid dataset or select a predefined option.")
    } else {
      return(NULL)
    }
  })
  
  # Render data table for the selected dataset
  output$lut_header <- DT::renderDataTable({
    dataset <- processed_data()
    if (!is.null(dataset) && ncol(dataset) > 0) {
      DT::datatable(
        dataset,
        options = list(
          scrollX = TRUE,
          columnDefs = list(list(visible = TRUE, targets = 0:4)),
          pageLength = 5
        ),
        selection = 'none'
      )
    } else {
      data.frame(Headers = c("No headers found in the dataset."))
    }
  })
  
  ## Get statiistical   ---------------------------------------
  
  # Reactive to calculate summary statistics for the selected bands
  band_stats <- reactive({
    
    validate(need(processed_data(), "No dataset available."))
    dataset <- processed_data()
    
    # Select only columns that start with 'B' (e.g., B1, B2, etc.)
    
    # Step 1: Select the bands (including B8A)
    band_cols <- dataset %>%
      dplyr::select(matches("^B[0-9]+[A-Z]?$"))
    
    
    if (input$plotting_indices == T ){
      # Show the modal window
      show_modal_spinner()
      #print(indices)
      showNotification("Generating spectral indices for model calculations...", type = "message")
      
      # Show a progress bar while processing
      withProgress(message = 'Getting spectral indices ...', {
        # Update progress
        total_steps <- 2
        incProgress(1/total_steps, detail = "")
        indices <-getIndicesSE2.ML(df=as.data.frame(band_cols), sensor = "Sentinel-2a", df.data = NULL, fast.process =T)
        
        # Final step
        incProgress(1/total_steps, detail = "Finalizing ...")
      })
      #print(indices)
      showNotification("Spectral indices generation completed.", type = "message")
      
      
      remove_modal_spinner()
      # After processing, update to the Indices tab
      updateTabsetPanel(session, "tabs", selected = "Indices")
      
      # Scroll to the top of the page using JavaScript
      session$sendCustomMessage(type = 'scrollToTop', message = list())
      
      if (is.null(indices) || nrow(indices) == 0) {
        showNotification("Error calculating indices. Please check your input data.", type = "error")
        return(NULL)
      }
      
      dataset_indices <- cbind(dataset, indices)
      showNotification("Indices calculated successfully.", type = "message")
    } else {
      showNotification("Please check the 'Calculate indices' box to proceed.", type = "warning")
      
    }
    
    
    # Calculate average, median, and percentiles for each band
    
    stats.SE2 <- band_cols %>%
      reframe(
        average = apply(band_cols, 2, mean, na.rm = TRUE),
        median = apply(band_cols, 2, median, na.rm = TRUE),
        percentile.25 = apply(band_cols, 2, quantile, probs = 0.25, na.rm = TRUE),
        percentile.50 = apply(band_cols, 2, quantile, probs = 0.50, na.rm = TRUE),
        percentile.75 = apply(band_cols, 2, quantile, probs = 0.75, na.rm = TRUE)
      ) %>%
      dplyr::mutate(Band = names(band_cols)) %>%
      dplyr::select(Band, everything())   %>%
      dplyr::mutate(Band = names(band_cols)) %>%
      dplyr::select(Band, everything()) %>%
      dplyr::filter(Band != "B8A" )
    # Sort the Band column as a factor
    stats.SE2 <- stats.SE2 %>%
      mutate(Band = factor(Band, levels = paste('B', sort(as.numeric(gsub('B', '', Band))), sep = "")))
    
    
    print(unique(stats.SE2$Band))
    showNotification(" Testing scatter-plot done successfully.", type = "warning")
    # Reshape data to long format for easier plotting
    stats_long <- stats.SE2 %>%
      pivot_longer(cols = -Band, names_to = "Statistic", values_to = "Value")
    
    if (input$plotting_indices == T ){
      
      return(list(stats_long = stats_long, dataset_indices = dataset_indices))
    } else {
      return(list(stats_long = stats_long))
    }
    
  })
  ## Get plot_checking   ---------------------------------------
  
  # Render the plot in UI
  output$plot_checking <- renderPlot({
    stats_long <- band_stats()$stats_long
    # Color if I need
    color.selected <- c('forestgreen','red','darkgoldenrod2','dodgerblue3','darkseagreen3','indianred2')
    color.i = color.selected[4]
    
    # Generate the plot using ggplot2
    plot_ <- ggplot(stats_long, aes(x = Band, y = Value, color=Statistic, group = Statistic)) +
      geom_line() +  # Lines connecting points
      geom_point() +  # Points at each value
      
      labs(title = "Average, Median, and Percentiles by Band",
           x = "Band",
           y = "Value") +
      theme_bw() +
      labs(title = "", x = "",y = "Reflectance") +
      
      theme(
        text = element_text(size = 14, face='bold'),  # Increase the text size
        axis.title = element_text(size = 16, face = "bold"),  # Make axis titles bold
        axis.text = element_text(face = "bold"),  # Make axis numbers bold
        axis.text.x = element_text(angle = 45, hjust = 1,face='bold'),
        plot.title = element_text(face = "bold"),  # Make plot title bold
        plot.subtitle = element_text(size = 14,face='bold')  # Adjust subtitle size
      )
    
    plot_
    
    
  })
  # In your server function
  output$instructions_ui <- renderUI({
    if (input$plotting_indices && !is.null(band_stats()))  {  # Check if indices have been estimated
      tagList(
        br(),
        h4("Instructions:"),
        p("Select the variables for the X and Y axes below to explore their relationship."),
        p("Note₁: Use the dropdowns to select the variables you wish to plot."),
        p(em("Note₂: Ensure that the indices have been calculated.")),
        p(HTML("If the option is not available, please check the corresponding <strong>Calculate Indices</strong> box."))
      )
    }
  })
  
  ## Get scatter_plot   ---------------------------------------
  #Render scatter plot for exploring relationships between indices
  output$scatter_plot <- renderPlot({
    dataset_indices <- band_stats()$dataset_indices
    req(dataset_indices)  # Ensure dataset_indices is available
    
    x_var <- input$x_axis  # User-selected variable for the x-axis
    y_var <- input$y_axis  # User-selected variable for the y-axis
    
    # Create scatter plot if both x and y variables are selected
    if (!is.null(x_var) && !is.null(y_var)) {
      ggplot(dataset_indices, aes_string(x = x_var, y = y_var)) +
        geom_point(color = 'blue') +
        labs(title = paste("Scatter-plot of", y_var, "vs", x_var),
             x = x_var,
             y = y_var) +
        theme_bw() +
        theme(
          text = element_text(size = 14, face = 'bold'),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(face = "bold")
        )
    } else {
      ggplot() + labs(title = "Select variables to plot")  # Show empty plot with a message
    }
  })
  
  # UI for selecting variables for the scatter plot
  output$select_axes_ui <- renderUI({
    dataset_indices <- band_stats()$dataset_indices
    req(dataset_indices)  # Ensure dataset_indices is available
    
    # Create selectInput for X-axis and Y-axis
    tagList(
      selectInput("x_axis", "Select X-axis variable:", choices = names(dataset_indices), selected = NULL),
      selectInput("y_axis", "Select Y-axis variable:", choices = names(dataset_indices), selected = 'NDVI')
    )
  })
  # Set the file size limit for uploads
  options(shiny.maxRequestSize = 10 * 1024^2)  # 10 MB
  
  # Reactive output to indicate if a dataset has been uploaded
  output$hasUploaded <- reactive({
    return(!is.null(user_lut_data()))
  })
  
  outputOptions(output, "hasUploaded", suspendWhenHidden = FALSE)
  
  
  # Render links based on R's Notebook selection (notebook_r)
  output$lut_selection <- renderUI({
    if (input$lut_db == "upload") {
      tagList(
        # Output for the plot
        
        br(),
        # Description for the DataTable
        p('The selected dataset contains relevant information and visualizations.'),
        p('Please find the details below.'),
      )
    } else if (input$lut_db == "PROSAIL") {
      tagList(
        
        br(),
        # Description for the DataTable
        p('The selected PROSAIl dataset contains relevant information and visualizations.'),
        p('Please find the details below.'),
        
        
      )
    } else if (input$lut_db == "INFORM") {
      tagList(
        br(),
        
        # Description for the DataTable
        p('The selected INFORM dataset contains relevant information and visualizations.'),
        p('Please find the details below.'),
        
      )
    } else {
      NULL
    }
    
  })
  
  # Listen for changes in the LUT database selection input
  observeEvent(input$lut_db, {
    # Switch to the 'LUT' tab when the input changes
    updateTabsetPanel(session, inputId = "tabs", selected = "LUTs")
  })
  
  
  ## Get depVar   ---------------------------------------
  # Function to get the input bands based on the dependent variable
  get_input_bands <- function(depVar) {
    if (depVar %in% c('Cab')) {
      return(c('B2','B3','B4','B5','B6','B7','B8','CR.red.nir.1',
               'CR.red.nir.6','NDVI','TCARI_OSAVI','RedEg1','CIre','PSSRa','MCARI','OSAVI','GM1','Datt1',
               'NDRE','IRECI','CIgreen', 'CR.Brown'))
    } else if (depVar == 'LAI') {
      return(c('B2','B3','B4','B5','B6','B7','B8','B11','B12',
               'CR.red.nir.6','NDVI','TCARI_OSAVI','RedEg1','CIre','PSSRa','MCARI','OSAVI','GM1',
               'NDRE','IRECI','CIgreen','CR.SWIR'))
    } else if (depVar == 'leaf water content') {
      return(c('B8A','B11','B12','WDRVI','MNDVI','NDWI','NDWI2','SBI','WET','GVI','CR.SWIR'))
    } else if (depVar == 'Cbrown') {
      return(c('B2','B3','B4','B5','B6','B7','B8','EVI','TCARI',
               'CR.red.nir.6','NDVIv','kNDVI','TCARI_OSAVI','RedEg1','CIre','PSSRa','MCARI','OSAVI','GM1','Datt1',
               'NDRE','IRECI','CIgreen', 'CR.Brown'))
    }  else if (depVar == 'EWT') {
      return(c('B8','B11','B12','WDRVI','MNDVI','NDWI','NDWI2','SBI','WET','GVI','CR.SWIR'))
    }
    
  }
  
  
  # Function to get the input bands based on the dependent variable
  get_selected_input <- function(inputs_toML) {
    if (inputs_toML %in% c('indices')) {
      return(c('CR.red.nir.1',
               'CR.red.nir.6','NDVI','TCARI_OSAVI','RedEg1','CIre','PSSRa','MCARI','OSAVI','GM1','Datt1',
               'NDRE','IRECI','CIgreen', 'CR.Brown'))
      
    }  else if (inputs_toML == 'visible_indices') {
      return(c('CR.red.nir.1','CR.red.nir.6','NDVI','TCARI_OSAVI','RedEg1','CIre','PSSRa','MCARI','OSAVI','GM1','Datt1',
               'NDRE','IRECI','CIgreen', 'CR.Brown'))
      
    } else if (inputs_toML == 'rededge_indices') {
      return(c('CR.red.nir.1', 'CR.red.nir.6','RedEg1','CIre','CIgreen', 'CR.Brown'))
      
    } else if (inputs_toML == 'swir_indices') {
      return(c('WDRVI','MNDVI','NDWI','NDWI2','SBI','WET','GVI','CR.SWIR'))
      
    } else if (inputs_toML == 'reflectance') {
      return(c('B2','B3','B4','B5','B6','B7','B8','B11','B12'))
      
    } else if (inputs_toML == 'vnir_reflectance') {
      return(c('B2','B3','B4','B5','B6','B7','B8'))
      
    } else if (inputs_toML == 'swir_reflectance') {
      return(c('B8','B11','B12'))
      
    } else if (inputs_toML == 'indices_reflectance') {
      return(c('B2','B3','B4','B5','B6','B7','B8','B11','B12','CR.red.nir.1',
               'CR.red.nir.6','NDVI','TCARI_OSAVI','RedEg1','CIre','PSSRa','MCARI','OSAVI','GM1','Datt1',
               'NDRE','IRECI','CIgreen', 'CR.Brown','WDRVI','MNDVI','NDWI','NDWI2','SBI','WET','GVI','CR.SWIR'))
      
    } else if (inputs_toML == 'indices_reflectance_vnir') {
      return(c('B2','B3','B4','B5','B6','B7','B8','CR.red.nir.1',
               'CR.red.nir.6','NDVI','TCARI_OSAVI','RedEg1','CIre','PSSRa','MCARI','OSAVI','GM1','Datt1',
               'NDRE','IRECI','CIgreen', 'CR.Brown'))
      
    } else if (inputs_toML == 'indices_reflectance_swir') {
      return(c('B8','B11','B12','WDRVI','MNDVI','NDWI','NDWI2','SBI','WET','GVI','CR.SWIR'))
    }
    
    
  }
  
  
  
  observeEvent(input$train_model, {
    
    # Eliminar variables globales antes de entrenar un nuevo modelo
    #rm(list = ls(globalenv()), envir = globalenv())
    #gc()  # Forzar la recolección de basura
    # Extract input values
    
    # Check if the checkbox is not selected
    if (!input$plotting_indices) {
      showModal(modalDialog(
        title = "Indices Calculation Required",
        "Please check 'Calculate and plot vegetation indices' to proceed with the calculation of indices.",
        easyClose = TRUE,
        footer = NULL
      ))
      
      # Exit the event if the checkbox is not checked
      return()
    }
    
    
    # After training is done, update to the Predictions tab
    updateTabsetPanel(session, "tabs", selected = "Predictions")
    
    # Scroll to the top of the page using JavaScript
    session$sendCustomMessage(type = 'scrollToTop', message = list())
    
    model <- input$models_
    depVar <- input$depVar
    inputML <- input$inputs_toML
    
    # Define the data frame and dependent variable
    dataset <- band_stats()$dataset_indices
    
    # Check the LUT database selection
    if (is.null(dataset)) {
      showModal(modalDialog(
        title = "Dataset Selection Required",
        "Please select a valid LUT database or upload a new dataset before proceeding.",
        easyClose = TRUE,
        footer = NULL
      ))
      
      # Exit the event if no valid LUT database is selected
      return()
    }
    
    nrows <- nrow(dataset)
    print(nrows)
    
    # Check if the number of rows is less than 30
    if (nrows < 100) {
      # Check if the input sample size is not equal to 30
      if (input$p_samplesML < 80) {
        showModal(modalDialog(
          title = "Insufficient Data for training a model",
          paste("The dataset provided has only", nrows, "row(s).",
                "Please provide at least 80% of the data zxfor training the model."),
          
          easyClose = TRUE,
          footer = NULL
        ))
        return()
      }
      
    }
    
    if (input$depVar == 'EWT' | input$depVar == 'Prot'  | input$depVar == 'CBC' ){
      dataset[[depVar]] <- dataset[[depVar]] * 1000
    }
    #inputs_bands <- get_input_bands(depVar)
    inputs_bands <- get_selected_input(inputML)
    print(inputs_bands)
    # Generate a random seed based on the current time
    set.seed(as.numeric(Sys.time()) %% 10000)  # Limits seed to 4 digits
    
    # Define the directory for model outputs
    output_dir <- file.path('www', 'models.ML')
    
    # Remove the existing directory if it exists
    if (dir.exists(output_dir)) {
      unlink(output_dir, recursive = TRUE)
    }
    # Create a temporary directory
    # Define the directory for model outputs
    temp_dir <- file.path(output_dir, format(Sys.time(), "%Y%m%d_%H%M%S"),'/')
    dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Generate filename with date and time
    paste0("model_and_stats_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tar")
    # Convert percentage to proportion
    n.prop <- input$p_samplesML / 100
    # Calculate the number of samples to reduce
    n.samples.reduced <- ceiling(n.prop * nrow(dataset))
    # Sample rows
    rows.r <- sample(nrow(dataset), n.samples.reduced)
    
    # Set a minimum row requirement
    min_rows <- 5  # You can adjust this based on your needs
    if (length(rows.r) <= min_rows) {
      # Show error notification and adjust the proportion to 10%
      showNotification(
        paste("Error: This dataset has very few simulations. Please consider providing a different LUT dataset."),
        type = "error",
        duration = 5
      )
      
      # Adjust proportion to 10% and recalculate
      n.prop <- 10 / 100  # Adjusted to 10% (was 20% in the original code)
      n.samples.reduced <- ceiling(n.prop * nrow(dataset))
      rows.r <- sample(nrow(dataset), n.samples.reduced)
      print(length(rows.r))
      # Show warning about the change in proportion
      showNotification(
        paste("Warning: The number of data samples for training and testing has been increased to 10%."),
        type = "warning",
        duration = 5
      )
    } else {
      # Convert percentage to proportion
      n.prop <- input$p_samplesML / 100
      # Calculate the number of samples to reduce
      n.samples.reduced <- ceiling(n.prop * nrow(dataset))
      # Sample rows
      rows.r <- sample(nrow(dataset), n.samples.reduced)
      print(length(rows.r))
    }

    show_modal_spinner()
    
    if (input$models_ =='CNN' | input$models_ == 'Hidden_layers'){
      
      if (input$models_ %in% c('CNN', 'Hidden_layers')) {
        model_ <- if (input$models_ == 'Hidden_layers') 'Hidden-layers' else input$models_
      }
      print(model_)
      
      showNotification("Loading the Python environment for Keras requires time the first time. Please be patient...", type = "warning") # Use withProgress to show a progress bar during the model training process
      withProgress(message = 'Training model', value = 0, {
        
        method.preProcess <- input$method.preProcess
        optimizer <- input$optimizer
        n_layers <- input$n_layers
        n_neurons <- input$n_neurons
        batch_size <- input$batch_size
        n_epochs <- input$n_epochs
        # Call the getMLmodel.withRetrain function
        models.time <- getMLmodel.withRetrain_app(
          dataset = dataset[rows.r, c(depVar, inputs_bands)],
          depVar = depVar,
          model = model_,
          optimizer = optimizer,
          n.times = 1,
          n.neurons = n_neurons,
          n.layers = n_layers,
          batch.size = batch_size,
          n.epochs = n_epochs,
          save.model = T,
          path.model = temp_dir,
          prop.split = c(0.8, 0.2),
          data.trans = 'preProcess',
          method.preProcess = method.preProcess,
          depVar.trans = FALSE,
          session=session )
      }) ## end training the model
      
      
      # Render the training history table
      #output$training_history <- DT::renderDT({
      # DT::datatable(history_df)
      #})
      output$status <- renderText("training model completed sucessfully.")
      
    } else {
      
      ## Best RF; SVM (no importance); NN; GB, xGB;
      algorithms<-c('RF','SVM','NN','GB', 'xGB','BRNN','Ensemble')
      algorithm.i  <- input$models_
      
      showNotification("Trainning the selected ML model...", type = "message")
      
      # Show a progress bar while processing
      withProgress(message = 'Trainning the ML model ...', {
        # Update progress
        total_steps <- 2
        incProgress(1/total_steps, detail = "")
        models.time<- get.inversion(data=dataset[rows.r, c(depVar, inputs_bands)], depVar=depVar, inputs=inputs_bands,n.cores=2,
                                    n.samples=n.samples.reduced,algorithm=algorithm.i,method.resampling='repeatedcv',
                                    seed=set.seed(as.numeric(Sys.time()) %% 10000) , save.model = T, save.path = temp_dir)
        
        # Final step
        incProgress(1/total_steps, detail = "Finalizing ...")
      })
      
      showNotification("Trainning selected ML model completed successfully.", type = "message")
      
      # print(models.time)
      
    }
    
    
    # Print the model summary to the viewer
    # output$model_summary <- renderPrint({
    #  summary(models.time[['models.keras']])
    #})
    
    # Define the download handler
    output$downloadData_ML <- downloadHandler(
      filename = function() {
        paste0("model_and_stats_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tar")
      },
      content = function(file) {
        
        current_date_dir <- temp_dir # Updated to current date
        # List all files in the current date directory
        files_to_compress <- list.files(path = current_date_dir, full.names = TRUE)
        
        # Print files to console (for debugging)
        print(files_to_compress)
        
        # Check if there are files to compress
        if (length(files_to_compress) > 0) {
          
          #Define the name of the output tar file
          tar_file_name <- paste0("www/downloads/model_and_stats_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tar")
          
          # Create the tar file, compressing all listed files
          tar(tar_file_name, files = files_to_compress)
          
          # Output a message indicating successful compression
          message(paste("Compressed files into:", tar_file_name))
          
          # Remove all files and folders in the www/models.ML/ directory
          ML_outputs <- 'www/models.ML/'
          all_folders_to_rm <- list.dirs(ML_outputs, full.names = TRUE, recursive = FALSE)  # List only the immediate subdirectories
          unlink(all_folders_to_rm, recursive = TRUE, force = TRUE)  # Remove only the folders
          message("Removed all folders in the www/models.ML/ directory.")
          
          # Additional check for the downloads directory
          downloads_dir <- 'www/downloads/'
          # if (length(list.files(downloads_dir)) > 2) {
          #  unlink(file.path(downloads_dir, "*"), recursive = TRUE, force = TRUE)  # Remove all files in downloads
          # message("Removed all files in the downloads directory.")
          # }
          
        } else {
          message("No files to compress in the specified directory.")
        }
        
        # Move the tar file to the specified output file location
        file.copy(tar_file_name, file)
      }
    )
    
    # Create reactive values to manage readiness state
    output_ready <- reactiveValues(plot = FALSE, table = FALSE)
    
    # Generate a ggplot (for demonstration)
    output$predictionPlot <- renderPlot({
      
      if (input$models_ =='CNN' | input$models_ == 'Hidden_layers'){
        plot_val <- models.time[['plot.val']]
      } else {
        plot_val <- models.time[['plot']]
      }
      print(plot_val)
      
      output_ready$plot <- TRUE
    })
    
    # Show notification after model training is complete
    observe({
      if (input$models_ =='CNN' | input$models_ == 'Hidden_layers'){
        if (!is.null(models.time[['plot.val']]) && !is.null(models.time[['stats']])) {
          showNotification("training model completed sucessfully.", type = "message")
        }
      } else {
        if (!is.null(models.time[['plot']]) && !is.null(models.time[['statistics']])) {
          # showNotification("Model training complete.", type = "message")
        }
      }
    })
    
    
    # Render table for the dataset in models.time[['stats']]
    output$statsTable_ML <- DT::renderDataTable({
      # Access the dataset
      if (input$models_ =='CNN' | input$models_ == 'Hidden_layers'){
        stats_df <- models.time[['stats']]
      } else {
        stats_df <- models.time[['statistics']]
      }
      
      
      # Render the dataset as a table
      DT::datatable(stats_df, options = list(pageLength = 10, autoWidth = TRUE))
      
      # output_ready$table <- TRUE
    })
    
    
    remove_modal_spinner()
    
    
    # Reactive block to generate code based on user selections
    model_code <- reactive({
      
      # If the selected model is CNN or Hidden Layers Neural Network
      if (input$models_ == 'CNN' | input$models_ == 'Hidden_layers') {
        
        # Define the model name
        model_ <- if (input$models_ == 'Hidden_layers') 'Hidden-layers' else input$models_
        
        # Generate code for CNN/Hidden Layers Neural Network model
        code <- paste0(
          "# Define the selected model and key hyperparameters\n",
          "model_ <- '", model_, "'\n",  # Model type (CNN or Hidden Layers)
          "method.preProcess <- '", input$method.preProcess, "'\n",  # Data pre-processing method
          "optimizer <- '", input$optimizer, "'\n",  # Optimizer type for model training
          "n_layers <- ", input$n_layers, "\n",  # Number of layers in the neural network
          "n_neurons <- ", input$n_neurons, "\n",  # Number of neurons per layer
          "batch_size <- ", input$batch_size, "\n",  # Batch size for training
          "n_epochs <- ", input$n_epochs, "\n",  # Number of epochs for training
          "n.prop <- ", input$p_samplesML_keras, " / 100\n",  # Proportion of data to use for training
          # Properly format inputs_bands as a vector
          "inputs_bands <- c(",paste("'",input$inputs_bands, "'",sep='', collapse = ", "), ")\n",  # Input bands selected by the user
          
          "# Call to function to train and evaluate the model\n",
          "models.time <- getMLmodel.withRetrain(\n",
          "  dataset = dataset[, c(",input$depVar,", inputs_bands)],\n",  # Use dependent variable and selected input bands
          "  depVar = depVar,\n",  # Dependent variable (trait to estimate)
          "  model = model_,\n",  # Model type (CNN/Hidden Layers)
          
          "  optimizer = optimizer,\n",  # Optimizer for training
          "  n.times = 1,\n",  # Number of times to repeat training
          "  n.neurons = n_neurons,\n",  # Number of neurons per layer
          "  n.layers = n_layers,\n",  # Number of layers
          "  batch.size = batch_size,\n",  # Batch size
          "  n.epochs = n_epochs,\n",  # Number of epochs
          "  save.model = TRUE,\n",  # Option to save the trained model
          "  path.model = temp_dir,\n",  # Path to save the model
          "  prop.split = c(0.8, 0.2),\n",  # Split ratio for training and testing data (80-20%)
          "  data.trans = 'preProcess',\n",  # Data transformation method
          "  method.preProcess = method.preProcess,\n",  # Pre-processing method
          "  depVar.trans = FALSE\n)"  # No transformation on dependent variable
          
        )
        
      } else {  # For other models (e.g., Random Forest, SVM, XGBoost, etc.)
        
        # Generate code for other machine learning models
        code <- paste0(
          "# Define the selected algorithm and key settings\n",
          "algorithm.i <- '", input$models_, "'\n",  # Selected algorithm (RF, SVM, XGBoost, etc.)
          "n.prop <- ", input$p_samplesML, " / 100\n",  # Proportion of data to use for training
          "n.samples.reduced <- ceiling(n.prop * nrow(dataset))\n",  # Calculate the number of samples for training
          # Properly format inputs_bands as a vector
          "inputs_bands <- c(", paste("'", get_selected_input(input$inputs_toML), "'", sep='',collapse = ", "), ")\n",  # Input bands selected by the user
          
          "# Call to function to train and evaluate the model\n",
          "models.time <- get.inversion(\n",
          "  data = dataset[, c(",input$depVar,", inputs_bands)],\n",  # Use dependent variable and selected input bands
          "  depVar = depVar,\n",  # Dependent variable (trait to estimate)
          "  inputs = inputs_bands,\n",  # Input bands or features selected by the user
          "  n.cores = 2,\n",  # Number of CPU cores to use for parallel processing
          "  n.samples = n.samples.reduced,\n",  # Number of samples for training
          "  algorithm = algorithm.i,\n",  # Algorithm type
          "  method.resampling = 'repeatedcv',\n",  # Resampling method (repeated cross-validation)
          "  seed = set.seed(as.numeric(Sys.time()) %% 10000),\n",  # Set a random seed for reproducibility
          "  save.model = TRUE,\n",  # Option to save the trained model
          "  save.path = temp_dir\n)"  # Path to save the model
          
        )
      }
      
      # Return the generated code as output
      return(code)
    })
    
    # Render the R code in the UI
    output$model_code <- renderText({
      model_code()
    })
    
    output$download_code <- downloadHandler(
      filename = function() { paste0(input$models_, "_model_code.R") },
      content = function(file) {
        writeLines(model_code(), file)
      }
    )
    
    # Check the ML is done, in this context downloadbutton apparece
    output$plotReady_tableReady <- renderUI({
      # Ensure both plot and table are ready
      req(output_ready$plot)  # Ensure plot data is available
      
      tagList(
        div(
          style = "display: inline-block; margin-right: 10px;",  # Style for first button
          downloadButton("downloadData_ML", "Download Model & Stats")
        ),
        
        div(
          style = "display: inline-block;",  # Style for second button
          downloadButton("download_code", "Download R Code")
        ),
        br(), br(),
        h4("Dowloading also the R Code used in your ML parametrization:"),
        verbatimTextOutput("model_code"),  # Display the generated R code
        # Properly format inputs_bands as a vector
        #inputs_bands <- paste("'", get_selected_input(input$inputs_toML), "'", sep = '', collapse = ", "),
        
        
        
        p(em('Note₁: The download includes the trained model and associated statistical data.')),
        p(em('Note₂: The R download button provides the main script for model training.')),
        br(),
        p(("The inputs used in the selected model are: ")),  # Subscript note for input bands
        DT::dataTableOutput("inputsBandsTable"),  # Placeholder for the inputs_bands DT table  # Placeholder for the inputs_bands table
        
        
      )
    })
    
    # Render the inputs_bands table using DT
    output$inputsBandsTable <- DT::renderDataTable({
      # Create a data frame with inputs_bands for display
      inputs_bands_vector <- get_selected_input(input$inputs_toML)  # Get selected input bands
      inputs_df <- data.frame(Inputs = inputs_bands_vector)  # Create a data frame
      # Render the dataset as a table
      DT::datatable(inputs_df, options = list(pageLength = 10, autoWidth = TRUE))
      
    })
    
  }) ## ed observer
  
}

# Create Shiny object
shinyApp(ui = ui, server = server)

