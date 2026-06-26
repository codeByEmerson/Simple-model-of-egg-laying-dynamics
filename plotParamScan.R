# handles plotting of output from gridscanner.R

#  RColorBrewer code derived from R-graph-gallery
# https://r-graph-gallery.com/40-rcolorbrewer-get-a-longer-palette.html
library(RColorBrewer)

plotpointsize=3#1.5#3.4

# to avoid calculting this manually in multiple places
# returns the maximum allowable value when rescaling for ease of comparison 
getErrorCap <- function(errorList, errorType="other") {    
    if (errorType == "mean") {
        return(5)
    } else if (errorType =="hist") {
        return(100)
    } else if (errorType == "megamix") { # mean+hist+slope+variance error combined
        #print(errorList)
        #print(median(errorList))
        return(4) # since errors are rescaled to errorclass/median(errorclass), 1 is a good cutoff per errorclass we are combining, since we combine 4, 4 is a good net cutoff. not perfect but reasonable.
        
        #return(median(errorList))
    } else {
        # use median to help avoid tail-skew in the mean, gives nice results
        return(median(errorList) * 2)
    }
}

# maps error values to colors for plotting
colormap_meanError <- function(errorList, errorType="other") {
    errorCap <- getErrorCap(errorList, errorType)
    scaleError <- errorList/errorCap

    # orange to blue
    palet <- brewer.pal(4, "PuOr")  

    # paletrange = # of colors total
    paletrange<-25

    # Add more colors to this palette :
    palet <- colorRampPalette(palet)(paletrange)
    # flip it so that low is blue, high is orange (orange = bad more intuitive)
    palet <- rev(palet)
    
    # add black as over-saturation indicator
    palet <- c(palet, "#000000")
    
    # map error values to palet values
    # 0-1 error value * palet length -> palet index for corresponding color
    # use CEILING and not ROUND to avoid index(0) errors !!!
    paletScaledError <- ceiling(scaleError * length(palet))
    
    # map oversaturated values (>errorCap) to final index (black)
    # without this step, oversaturated's are white
    paletScaledError[paletScaledError>length(palet)]<-length(palet)
    
    # turn the vector of errors into a vector of colors corresponding to those errors
    colorified <- palet[paletScaledError]
    
    return(colorified)
}

# deprecated.
make3Dplot <- function(paramX, paramY, paramZ, name="3dplot", errTrim=TRUE) {
    library(plotly)

    # for 3d plot, you need 1 meanError for every As/Ad combination
    # so length(meanError) = length(As)*length(Ad)
    # to do this, we build and populate a matrix
    
    # to avoid hugely bad things skewing perspective
    if (errTrim) {
        # less than b/c we are in negative space here
        paramZ[paramZ<median(paramZ)*2]<-median(paramZ)*2
    }
    
    # prep dataframe
    xyzdat <- data.frame(
        'x' = paramX,
        'y' = paramY,
        'z' = paramZ
    )
    
    # convert to matrix 
    xyzmat<-xtabs(z ~ x + y, data = xyzdat)
    
    # assign rownames so the plot has correct axis labels
    cx <- as.numeric(colnames(xyzmat))
    cy <- as.numeric(rownames(xyzmat))



    q<- plot_ly(z = ~(xyzmat), x=cx, y=cy) %>% add_surface() 
    # derived from: https://community.plotly.com/t/defined-axes-with-limits/6126
    q<-layout(q, title = 'Manually Specified Labels',
                 scene = list(
                                xaxis = list(title = "Activator half-life (A1/2)"),   # Change x/y/z axis title
                                yaxis = list(title = "Activator release size (As)"),
                                zaxis = list(title = "Difference from experimental data")))
    


    prefixName <- paste("/htmlout/", name)
    postfixName <- paste(prefixName, ".html")

    library(htmlwidgets)
    saveWidget(q, file=paste0( getwd(), postfixName))
}

# handles plotting of output from gridscanner.R
plotParamScan <- function(uName, fitIhl=4, fitIs=7, fitAp=0.125, fitAd=6, fitAs=2, geno="unspecified") {

    csvfile<-paste(uName, ".csv")
    figName<-paste(uName, ".pdf")
    pdf(figName)
    trials_df <- read.csv(csvfile)

    trials_df_nona <- na.omit(trials_df)


    plot(density(trials_df_nona$meanInactiveLength),xlim=c(0,120))


    if (geno == "wt") {
        plot(trials_df_nona$param_Ihl,trials_df_nona$meanError, main="Model 2 (WT, complex), scan of Ihl values 0-6 versus\n RMSE of Mean Inactive phase",
        ylab="RMSE(Mean Inactive Phase Duration)",
        xlab="Ihl value (4 is value in paper)",
        ylim=c(0,10))
        
        abline(v=fitIhl, col="blue")

        
        plot(trials_df_nona$param_Ihl,trials_df_nona$histError, main="Model 2 (WT, complex), scan of Ihl values 0-6 versus\n RMSE of history dependence",
        ylab="RMSE(Mean Inactive Phase Duration)",
        xlab="Ihl value (4 is value in paper)",
        ylim=c(0,200))
        
        abline(v=fitIhl, col="blue")
        
        plot(trials_df_nona$param_Is,trials_df_nona$meanError, main="Model 2 (WT, complex), scan of Is values 0-14 versus\n RMSE of Mean Inactive phase",
        ylab="RMSE(Mean Inactive Phase Duration)",
        xlab="Is value (7 is value in paper)",
        ylim=c(0,10))
        
        abline(v=fitIs, col="blue")
        
        plot(trials_df_nona$param_Is,trials_df_nona$histError, main="Model 2 (WT, complex), scan of Is values 0-14 versus\n RMSE of history dependence",
        ylab="RMSE(Mean Inactive Phase Duration)",
        xlab="Is value (7 is value in paper)",
        ylim=c(0,200))
        
        abline(v=fitIs, col="blue")

        # plots distribution of errors and show color - error cutoff!!!    
        plot(density(trials_df_nona$meanError), main="RMSE(Mean inactive phase)")
        abline(v=getErrorCap(trials_df_nona$meanError, "mean"), col="black", lwd=3)
        
        plot(density(trials_df_nona$histError), main="RMSE(History Dependence)")
        abline(v=getErrorCap(trials_df_nona$histError, "hist"), col="black", lwd=3)

        plot(trials_df_nona$param_Ihl,trials_df_nona$param_Is, main="Model 2 (WT, complex), scan of Is values 0-14 versus\n Ihl values 0-6 versus mean inactive phase error (color)",
        col=colormap_meanError(trials_df_nona$meanError, "mean"),
        xlab="Ihl value (4 is value in paper)",
        ylab="Is value (7 is value in paper)",
        pch=15, cex=plotpointsize)
        
        abline(v=fitIhl, col="grey", lwd=3)
        abline(h=fitIs, col="grey", lwd=3)
        
        plot(trials_df_nona$param_Ihl,trials_df_nona$param_Is, main="Model 2 (WT, complex), scan of Is values 0-14 versus\n Ihl values 0-6 versus history error (color)",
        col=colormap_meanError(trials_df_nona$histError, "hist"),
        xlab="Ihl value (4 is value in paper)",
        ylab="Is value (7 is value in paper)",
        pch=15, cex=plotpointsize)
        
        abline(v=fitIhl, col="grey", lwd=3)
        abline(h=fitIs, col="grey", lwd=3)
        
        plot(density(trials_df_nona$histError/max(trials_df_nona$histError)))
        plot(density(trials_df_nona$meanError/max(trials_df_nona$meanError)))
        plot(density(trials_df_nona$meanError/max(trials_df_nona$meanError) + trials_df_nona$histError/max(trials_df_nona$histError)))
        
        
        plot(trials_df_nona$param_Ihl,trials_df_nona$param_Is, main="Model 2 (WT, complex), scan of Is values 0-14 versus\n Ihl values 0-6 versus inactive + history error (color)\n (not quite working)",
        col=colormap_meanError(
            (trials_df_nona$histError/max(trials_df_nona$histError) + trials_df_nona$meanError/max(trials_df_nona$meanError))/2 ),
        xlab="Ihl value (4 is value in paper)",
        ylab="Is value (7 is value in paper)",
        pch=15,
        cex=plotpointsize)#19
        
        abline(v=fitIhl, col="grey", lwd=3)
        abline(h=fitIs, col="grey", lwd=3)
        
        
        plot(trials_df_nona$param_Ihl,trials_df_nona$param_Is, main="Model 2 (WT, complex), scan of Is versus\n Ihl versus Error (mean inactive+history) )",
        col=colormap_meanError(trials_df_nona$meanError, "mean"),
        xlab="Ihl value (4 is value in paper)",
        ylab="Is value (7 is value in paper)",
        pch=15, cex=plotpointsize)
        
        points(trials_df_nona$param_Ihl,trials_df_nona$param_Is,
            col=colormap_meanError(trials_df_nona$histError, "hist"),
            pch=15, cex=plotpointsize*.588)
        
        abline(v=fitIhl, col="grey", lwd=3)
        abline(h=fitIs, col="grey", lwd=3)
        
        ######
        ######
        ######
        
        plot(trials_df_nona$param_Ihl,trials_df_nona$param_Is, main="Model 2 (WT, complex), scan of Is values versus\n Ihl versus slopeError",
        col=colormap_meanError(trials_df_nona$model30SlopeError),
        ylab="Is value",
        xlab="Ihl value",
        pch=15,
        cex=2)
        
        plot(trials_df_nona$param_Ihl,trials_df_nona$param_Is, main="Model 2 (WT, complex), scan of Is values versus\n Ihl versus histError+meanError and slopeError",
        col=colormap_meanError(trials_df_nona$meanError/median(trials_df_nona$meanError)+trials_df_nona$histError/median(trials_df_nona$histError)),
        ylab="Is value ",
        xlab="Ihl value ",
        pch=15, cex=plotpointsize)
        
        points(trials_df_nona$param_Ihl,trials_df_nona$param_Is,
        col=colormap_meanError(trials_df_nona$model30SlopeError),
        pch=15, cex=plotpointsize*.588)
        
        abline(v=fitIhl, col="grey", lwd=3)
        abline(h=fitIs, col="grey", lwd=3)
        
        plot(trials_df_nona$param_Ihl,trials_df_nona$param_Is, main="Model 2 (WT, complex), scan of Is values versus\n Ihl versus varError",
        col=colormap_meanError(trials_df_nona$varError),
        ylab="Is value",
        xlab="Ihl value",
        pch=15,
        cex=2)
        
        abline(v=fitIhl, col="grey", lwd=3)
        abline(h=fitIs, col="grey", lwd=3)
        
        
        plot(trials_df_nona$param_Ihl,trials_df_nona$param_Is, main="Model 2 (WT, complex), scan of Is values versus\n Ihl versus histError+meanError+slopeError+varError",
        col=colormap_meanError(trials_df_nona$meanError/median(trials_df_nona$meanError)+trials_df_nona$histError/median(trials_df_nona$histError) + trials_df_nona$model30SlopeError/median(trials_df_nona$model30SlopeError) + trials_df_nona$varError/median(trials_df_nona$varError), "megamix"),
        ylab="Is value ",
        xlab="Ihl value ",
        pch=15, cex=plotpointsize*0.75,
        ylim=c(0,14.7),
        xlim=c(0,17))
        
        plot(trials_df_nona$param_Ihl,trials_df_nona$param_Is, main="Model 2 (WT, complex), scan of Is values versus\n Ihl versus histError+meanError+slopeError+varError",
        col=colormap_meanError(trials_df_nona$meanError/median(trials_df_nona$meanError)+trials_df_nona$histError/median(trials_df_nona$histError) + trials_df_nona$model30SlopeError/median(trials_df_nona$model30SlopeError) + trials_df_nona$varError/median(trials_df_nona$varError), "megamix"),
        ylab="Is value ",
        xlab="Ihl value ",
        pch=15, cex=plotpointsize*0.75,
        ylim=c(0,14.7),
        xlim=c(0,17))
        
        abline(v=fitIhl, col="grey", lwd=3)
        abline(h=fitIs, col="grey", lwd=3)
        
        
        plot(trials_df_nona$param_Ad,trials_df_nona$param_As, main="Model 2 (WT, complex), scan of As values versus\n Ahl versus histError+meanError+slopeError+varError",
        col=colormap_meanError(trials_df_nona$meanError/median(trials_df_nona$meanError)+trials_df_nona$histError/median(trials_df_nona$histError) + trials_df_nona$model30SlopeError/median(trials_df_nona$model30SlopeError) + trials_df_nona$varError/median(trials_df_nona$varError), "megamix"),
        ylab="As value ",
        xlab="Ahl value ",
        pch=15, cex=plotpointsize*0.8,
        xlim=c(0,8.1),
        ylim=c(0,6.9))
        
        plot(trials_df_nona$param_Ad,trials_df_nona$param_As, main="Model 2 (WT, complex), scan of As values versus\n Ahl versus histError+meanError+slopeError+varError",
        col=colormap_meanError(trials_df_nona$meanError/median(trials_df_nona$meanError)+trials_df_nona$histError/median(trials_df_nona$histError) + trials_df_nona$model30SlopeError/median(trials_df_nona$model30SlopeError) + trials_df_nona$varError/median(trials_df_nona$varError), "megamix"),
        ylab="As value ",
        xlab="Ahl value ",
        pch=15, cex=plotpointsize*.8,
        xlim=c(0,8.1),
        ylim=c(0,6.9))
        
        abline(v=fitAd, col="grey", lwd=3)
        abline(h=fitAs, col="grey", lwd=3)
        
    }
        
    if (geno == "flp-11 nlp-7") {
            
        plot(trials_df_nona$param_Ihl,trials_df_nona$param_Is, main="Model 2 (FLP11NLP7, complex), scan of Is values versus\n Ihl versus histError+meanError+slopeError+varError",
        col=colormap_meanError(trials_df_nona$meanError/median(trials_df_nona$meanError)+trials_df_nona$histError/median(trials_df_nona$histError) + trials_df_nona$model30SlopeError/median(trials_df_nona$model30SlopeError) + trials_df_nona$varError/median(trials_df_nona$varError), "megamix"),
        ylab="Is value ",
        xlab="Ihl value ",
        pch=15, cex=plotpointsize*0.75,
        ylim=c(0,14.7),
        xlim=c(0,17))
        
        plot(trials_df_nona$param_Ihl,trials_df_nona$param_Is, main="Model 2 (FLP11NLP7, complex), scan of Is values versus\n Ihl versus histError+meanError+slopeError+varError",
        col=colormap_meanError(trials_df_nona$meanError/median(trials_df_nona$meanError)+trials_df_nona$histError/median(trials_df_nona$histError) + trials_df_nona$model30SlopeError/median(trials_df_nona$model30SlopeError) + trials_df_nona$varError/median(trials_df_nona$varError), "megamix"),
        ylab="Is value ",
        xlab="Ihl value ",
        pch=15, cex=plotpointsize*0.75,
        ylim=c(0,14.7),
        xlim=c(0,17))
        
        abline(v=fitIhl, col="grey", lwd=3)
        abline(h=fitIs, col="grey", lwd=3)
    }
    
     if (geno == "tdc-1") {

        plot(trials_df_nona$param_Ihl,trials_df_nona$param_Is, main="Model 2 (tdc-1, complex), scan of Is values versus\n Ihl versus histError+meanError+slopeError+varError",
        col=colormap_meanError(trials_df_nona$meanError/median(trials_df_nona$meanError)+trials_df_nona$histError/median(trials_df_nona$histError) + trials_df_nona$model30SlopeError/median(trials_df_nona$model30SlopeError) + trials_df_nona$varError/median(trials_df_nona$varError), "megamix"),
        ylab="Is value ",
        xlab="Ihl value ",
        pch=15, cex=plotpointsize*0.75,
        xlim=c(0,17))

        plot(trials_df_nona$param_Ihl,trials_df_nona$param_Is, main="Model 2 (tdc-1, complex), scan of Is values versus\n Ihl versus histError+meanError+slopeError+varError",
        col=colormap_meanError(trials_df_nona$meanError/median(trials_df_nona$meanError)+trials_df_nona$histError/median(trials_df_nona$histError) + trials_df_nona$model30SlopeError/median(trials_df_nona$model30SlopeError) + trials_df_nona$varError/median(trials_df_nona$varError), "megamix"),
        ylab="Is value ",
        xlab="Ihl value ",
        pch=15, cex=plotpointsize*0.75,
        xlim=c(0,17))

        abline(v=fitIhl, col="grey", lwd=3)
        abline(h=fitIs, col="grey", lwd=3)
    
    }
    
    if (geno == "tph-1") {
        
        plot(trials_df_nona$param_Ad,trials_df_nona$param_As, main="Model 2 (tph-1, complex), scan of As values versus\n Ahl versus histError+meanError+slopeError+varError",
        col=colormap_meanError(trials_df_nona$meanError/median(trials_df_nona$meanError)+trials_df_nona$histError/median(trials_df_nona$histError) + trials_df_nona$model30SlopeError/median(trials_df_nona$model30SlopeError) + trials_df_nona$varError/median(trials_df_nona$varError), "megamix"),
        ylab="As value ",
        xlab="Ahl value ",
        pch=15, cex=plotpointsize*.79,
        xlim=c(0,8))
        
        plot(trials_df_nona$param_Ad,trials_df_nona$param_As, main="Model 2 (tph-1, complex), scan of As values versus\n Ahl versus histError+meanError+slopeError+varError",
        col=colormap_meanError(trials_df_nona$meanError/median(trials_df_nona$meanError)+trials_df_nona$histError/median(trials_df_nona$histError) + trials_df_nona$model30SlopeError/median(trials_df_nona$model30SlopeError) + trials_df_nona$varError/median(trials_df_nona$varError), "megamix"),
        ylab="As value ",
        xlab="Ahl value ",
        pch=15, cex=plotpointsize*.79,
        xlim=c(0,8))
        
        
        abline(v=fitAd, col="grey", lwd=3)
        abline(h=fitAs, col="grey", lwd=3)
        
    }
    
    if (geno == "nlp-3") {
        
        plot(trials_df_nona$param_Ad,trials_df_nona$param_As, main="Model 2 (nlp-3, complex), scan of As values versus\n Ahl versus histError+meanError+slopeError+varError",
        col=colormap_meanError(trials_df_nona$meanError/median(trials_df_nona$meanError)+trials_df_nona$histError/median(trials_df_nona$histError) + trials_df_nona$model30SlopeError/median(trials_df_nona$model30SlopeError) + trials_df_nona$varError/median(trials_df_nona$varError), "megamix"),
        ylab="As value ",
        xlab="Ahl value ",
        pch=15, cex=plotpointsize*.8,
        xlim=c(0,8))
        
        plot(trials_df_nona$param_Ad,trials_df_nona$param_As, main="Model 2 (nlp-3, complex), scan of As values versus\n Ahl versus histError+meanError+slopeError+varError",
        col=colormap_meanError(trials_df_nona$meanError/median(trials_df_nona$meanError)+trials_df_nona$histError/median(trials_df_nona$histError) + trials_df_nona$model30SlopeError/median(trials_df_nona$model30SlopeError) + trials_df_nona$varError/median(trials_df_nona$varError), "megamix"),
        ylab="As value ",
        xlab="Ahl value ",
        pch=15, cex=plotpointsize*.8,
        xlim=c(0,8))
        
        abline(v=fitAd, col="grey", lwd=3)
        abline(h=fitAs, col="grey", lwd=3)
        
    }

}



