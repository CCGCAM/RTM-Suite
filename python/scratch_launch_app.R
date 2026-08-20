root <- "C:/Users/camin001/OneDrive - Wageningen University & Research/Workspace/0-RTM-Suite"
t0 <- Sys.time()
shiny::runApp(file.path(root, "AEO-Course/Apps/RTMs"), port = 8931, launch.browser = FALSE, host = "127.0.0.1")
