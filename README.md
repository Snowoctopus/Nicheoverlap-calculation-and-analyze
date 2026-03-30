# Niche overlap calculation
--> nicheoverlap.R is for the calculation of nicheoverlap

--> Running on Rstudio+R

--> Prepare data files: Use Matlab to produce the Depth+Phyto+AOA.txt file [writematrix(trab, 'trab.txt')]
    Data arrangement: 
    1st column-Depth; 2nd column-AOA abundance; 3th column-Phyto abundance 
    if use time-series data (e.g. k days)
    1st column-Depth; 2~k/2+1 column-AOA abundance; k/2+2~k+1 column Phyto abundance

-->Running command:
    source ('path/nicheoverlap.R') % load the function 
    trab <- read.table('path/trab.txt', sep=",", header=FALSE) % load the depth and abundance data in csv. file
    trab <- as.matrix(trab)
    NO <- nicheoverlap(trab)

--> Output:
    Nicheoverlap figure
    Nicheoverlap value recorded in the file defined by the 'sink' command
  

# FigurePlot_for-autotrophy-niche-analyze
These codes are designed for plotting figure in my submission paper about autotroph niche analyze
