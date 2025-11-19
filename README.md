# Sediment_Shape_Evolution
This repository contains all the code used to perform the analysis presented in our paper on modelling the evolution of sediment shape in the Sarzana Stream, submitted to the Journal of Geophysical Research: Earth Surface.
All colde are written in MATLAB and require the 'Signal Processing Toolbox' and the 'Statistics and Machine Learning Toolbox'.
## Repository structure
The file Main_Code.m controls the set of inputs, parameters and defines the output folders and files names.
In order to help interested users to reproduce the figures presented in the paper, individual folders are provided for each figure. If the input dataset is available, this is provided, together with a code to reproduce the presented image. If the input dataset is derived through an optimization procedure, the code to perform this optimization is provided within the REQUIRED_FUNCTIONS folder, and the code Main_Code.m controls its inputs and outputs locations and parameters. The definition of each parameter is provided in individual README files within each folder, together with the explanation on how to set-up the Main_Code file in order to obtain the desired figure.
