# Runs the model specified in model.R with various parameter choices

unknownGenoErr <- function () {
    print("ERROR ERROR ERROR ERROR ERROR ERROR")
    print("ERROR ERROR ERROR ERROR ERROR ERROR")
    print("ERROR ERROR ERROR ERROR ERROR ERROR")
    print("ERROR ERROR ERROR ERROR ERROR ERROR")
    print("ERROR ERROR ERROR ERROR ERROR ERROR")
    print("ERROR ERROR ERROR ERROR ERROR ERROR")
    print("ERROR ERROR ERROR ERROR ERROR ERROR")
    print("! genotype does not match scan type !")
}

gridscanner_ApAd <- function(OUTPUT_FILENAME, Ap_seq, Ad_seq, Is, Ihl, As, genotype, randN=6) {
    csvfilename <- paste(OUTPUT_FILENAME, ".csv")
    pdffilename <- paste(OUTPUT_FILENAME, ".pdf")

    # contains functions used in analysis of data
    source("analyzeParamScan.R")

    # contains the implementation of the model
    source("model.R")

    # contains code to run various parameter scans
    source("runParamScan.R")

    if (file.exists(csvfilename)) {
        lookupdf<-read.csv(csvfilename)
    } else {
        lookupdf<-data.frame()
    }

    # for storing all the trials
    trials<-vector("list", length(Ap_seq)*length(Ad_seq))

    # used for write indexing
    i<-1

    i.Ap<-1
    while (i.Ap<length(Ap_seq)) {

        i.Ad<-1
        while (i.Ad<length(Ad_seq)) {
            use_Is=Is
            use_As=As
            use_Ihl=Ihl
            use_Ad=Ad_seq[i.Ad]
            use_Ap=Ap_seq[i.Ap]
            use_capA=FALSE
            use_randIs=TRUE
            use_randAs=TRUE
            # == is unsafe b/c doubles acquire error, use isSame comparison instead
            # cant use all.equal b/c type issues (it wants numeric, these are double)
            alreadyDone<-any(isSame(lookupdf$param_Is,use_Is) & isSame(lookupdf$param_As,use_As) & isSame(lookupdf$param_Ihl,use_Ihl) & isSame(lookupdf$param_Ad,use_Ad) & isSame(lookupdf$param_Ap,use_Ap) & isSame(lookupdf$param_capA,use_capA) & isSame(lookupdf$param_randIs,use_randIs) & isSame(lookupdf$param_randAs,use_randAs) ) 
             

            if (alreadyDone == FALSE) {
                if (genotype == "tph-1") {
                    trials[[i]]<-runAndCalcTPH1_mod2(uIs=use_Is, uAs=use_As, uIhl=use_Ihl, uAd=use_Ad, uAp=use_Ap, ucapA=use_capA, urandIs=use_randIs, urandAs=use_randAs, urandN=randN)
                } else if (genotype == "nlp-3") {
                    trials[[i]]<-runAndCalcNLP3_mod2(uIs=use_Is, uAs=use_As, uIhl=use_Ihl, uAd=use_Ad, uAp=use_Ap, ucapA=use_capA, urandIs=use_randIs, urandAs=use_randAs, urandN=randN)
                } else
                {
                    unknownGenoErr()
                }
                
            }
            else {
                #print("skipping!")
            }
            i.Ad<-i.Ad+1
            i<-i+1
        }
        i.Ap<-i.Ap+1
    }

    trials <- trials[lengths(trials)!=0]

    # for rbindlist
    library('data.table')

    #turn our list of data.frames into one big data.frame
    trials_df <- rbindlist(trials)

    fixnames <- function(filename) {
        if (!file.exists(csvfilename)) {
            return(NA)
            # col.names=NA add a blank name to fix alignment
            # from: https://stackoverflow.com/questions/17890777/why-are-exported-row-names-always-off-by-1-in-r
        } else {
            return(FALSE) 
        }
    }

    write.table(trials_df, csvfilename, 
                sep = ",", col.names = fixnames(csvfilename), append = T)


    source('plotParamScan.R')
    if (genotype == "tph-1") {
        #Is=7, As=1.1, Ihl=4, Ad=0.15, Ap=0.0525, capA=FALSE, randIs=TRUE, randAs=TRUE, expDecay=TRUE, randN=3
        plotParamScan(uName=OUTPUT_FILENAME, 
                      fitIhl=4, fitIs=7, fitAp=0.0525, fitAd=4.5, fitAs=1.1, 
                      geno= "tph-1")
    } else if (genotype == "nlp-3") {
        # Is=7, As=0.5, Ihl=4, Ad=24, Ap=0.06, capA=FALSE, randIs=TRUE, randAs=TRUE, expDecay=TRUE, randN=1)
        plotParamScan(uName=OUTPUT_FILENAME, 
                      fitIhl=4, fitIs=7, fitAp=0.06, fitAd=1, fitAs=0.4, 
                      geno= "nlp-3")
    } else {
        unknownGenoErr()
    }
}



gridscanner_AsAd <- function(OUTPUT_FILENAME, As_seq, Ad_seq, Is, Ihl, Ap, genotype, randN=6) {
    csvfilename <- paste(OUTPUT_FILENAME, ".csv")
    pdffilename <- paste(OUTPUT_FILENAME, ".pdf")

    # contains functions used in analysis of data
    source("analyzeParamScan.R")

    # contains the implementation of the model
    source("model.R")

    # contains code to run various parameter scans
    source("runParamScan.R")

    if (file.exists(csvfilename)) {
        lookupdf<-read.csv(csvfilename)
    } else {
        lookupdf<-data.frame()
    }

    # for storing all the trials
    trials<-vector("list", length(As_seq)*length(Ad_seq))

    # used for write indexing
    i<-1

    i.As<-1
    while (i.As<length(As_seq)) {

        i.Ad<-1
        while (i.Ad<length(Ad_seq)) {
            use_Is=Is
            use_As=As_seq[i.As]
            use_Ihl=Ihl
            use_Ad=Ad_seq[i.Ad]
            use_Ap=Ap
            use_capA=FALSE
            use_randIs=TRUE
            use_randAs=TRUE
            # == is unsafe b/c doubles acquire error, use isSame comparison instead
            # cant use all.equal b/c type issues (it wants numeric, these are double)
            alreadyDone<-any(isSame(lookupdf$param_Is,use_Is) & isSame(lookupdf$param_As,use_As) & isSame(lookupdf$param_Ihl,use_Ihl) & isSame(lookupdf$param_Ad,use_Ad) & isSame(lookupdf$param_Ap,use_Ap) & isSame(lookupdf$param_capA,use_capA) & isSame(lookupdf$param_randIs,use_randIs) & isSame(lookupdf$param_randAs,use_randAs) ) 
             

            if (alreadyDone == FALSE) {
                if (genotype == "tph-1") {
                    trials[[i]]<-runAndCalcTPH1_mod2(uIs=use_Is, uAs=use_As, uIhl=use_Ihl, uAd=use_Ad, uAp=use_Ap, ucapA=use_capA, urandIs=use_randIs, urandAs=use_randAs, urandN=randN)
                } else if (genotype == "nlp-3") {
                    trials[[i]]<-runAndCalcNLP3_mod2(uIs=use_Is, uAs=use_As, uIhl=use_Ihl, uAd=use_Ad, uAp=use_Ap, ucapA=use_capA, urandIs=use_randIs, urandAs=use_randAs, urandN=randN)
                } else if (genotype == "wt") {
                    trials[[i]]<-runAndCalcWT_mod2(uIs=use_Is, uAs=use_As, uIhl=use_Ihl, uAd=use_Ad, uAp=use_Ap, ucapA=use_capA, urandIs=use_randIs, urandAs=use_randAs, urandN=randN)
                } else
                {
                    unknownGenoErr()
                }
                
            }
            else {
                #print("skipping!")
            }
            i.Ad<-i.Ad+1
            i<-i+1
        }
        i.As<-i.As+1
    }

    trials <- trials[lengths(trials)!=0]

    # for rbindlist
    library('data.table')

    #turn our list of data.frames into one big data.frame
    trials_df <- rbindlist(trials)

    fixnames <- function(filename) {
        if (!file.exists(csvfilename)) {
            return(NA)
            # col.names=NA add a blank name to fix alignment
            # from: https://stackoverflow.com/questions/17890777/why-are-exported-row-names-always-off-by-1-in-r
        } else {
            return(FALSE) 
        }
    }

    write.table(trials_df, csvfilename, 
                sep = ",", col.names = fixnames(csvfilename), append = T)


    source('plotParamScan.R')
    if (genotype == "tph-1") {
        #Is=7, As=1.1, Ihl=4, Ad=0.15, Ap=0.0525, capA=FALSE, randIs=TRUE, randAs=TRUE, expDecay=TRUE, randN=3
        plotParamScan(uName=OUTPUT_FILENAME, 
                      fitIhl=4, fitIs=7, fitAp=0.0525, fitAd=4.5, fitAs=1.2, 
                      geno= "tph-1")
    } else if (genotype == "nlp-3") {
        # Is=7, As=0.5, Ihl=4, Ad=24, Ap=0.06, capA=FALSE, randIs=TRUE, randAs=TRUE, expDecay=TRUE, randN=1)
        plotParamScan(uName=OUTPUT_FILENAME, 
                      fitIhl=4, fitIs=7, fitAp=0.06, fitAd=1, fitAs=0.4, 
                      geno= "nlp-3") 
    } else if (genotype == "wt") {
        # Is=7, As=0.5, Ihl=4, Ad=24, Ap=0.06, capA=FALSE, randIs=TRUE, randAs=TRUE, expDecay=TRUE, randN=1)
        plotParamScan(uName=OUTPUT_FILENAME, 
                      fitIhl=4, fitIs=7, fitAp=0.125, fitAd=1.75, fitAs=2, 
                      geno= "wt")
    } else {
        unknownGenoErr()
    }
}


gridscanner_IsIhl <- function(OUTPUT_FILENAME, As, Ad, Is_seq, Ihl_seq, Ap, genotype, randN=6) {
    csvfilename <- paste(OUTPUT_FILENAME, ".csv")
    pdffilename <- paste(OUTPUT_FILENAME, ".pdf")

    # contains functions used in analysis of data
    source("analyzeParamScan.R")

    # contains the implementation of the model
    source("model.R")

    # contains code to run various parameter scans
    source("runParamScan.R")

    if (file.exists(csvfilename)) {
        lookupdf<-read.csv(csvfilename)
    } else {
        lookupdf<-data.frame()
    }

    # for storing all the trials
    trials<-vector("list", length(Is_seq)*length(Ihl_seq))

    # used for write indexing
    i<-1

    i.Is<-1
    while (i.Is<length(Is_seq)) {

        i.Ihl<-1
        while (i.Ihl<length(Ihl_seq)) {

            #Is=7, As=2, Ihl=4, Ad=6, Ap=0.125, capA=FALSE, randIs=TRUE, randAs=TRUE, expDecay=TRUE

            use_Is=Is_seq[i.Is]
            use_As=As
            use_Ihl=Ihl_seq[i.Ihl]
            use_Ad=Ad
            use_Ap=Ap
            use_capA=FALSE
            use_randIs=TRUE
            use_randAs=TRUE
            # unsafe b/c doubles acquire error, use isSame comparison instead
            # cant use all.equal b/c type issues (it wants numeric, these are double)
         
            alreadyDone<-any(isSame(lookupdf$param_Is,use_Is) & isSame(lookupdf$param_As,use_As) & isSame(lookupdf$param_Ihl,use_Ihl) & isSame(lookupdf$param_Ad,use_Ad) & isSame(lookupdf$param_Ap,use_Ap) & isSame(lookupdf$param_capA,use_capA) & isSame(lookupdf$param_randIs,use_randIs) & isSame(lookupdf$param_randAs,use_randAs) ) 
             

            if (alreadyDone == FALSE) {
                if (genotype == "tdc-1") {
                    trials[[i]]<-runAndCalcTDC1_mod2(uIs=use_Is, uAs=use_As, uIhl=use_Ihl, uAd=use_Ad, uAp=use_Ap, ucapA=use_capA, urandIs=use_randIs, urandAs=use_randAs, urandN=randN)
                } else if (genotype == "flp-11 nlp-7") {
                    trials[[i]]<-runAndCalcFLP11NLP7_mod2(uIs=use_Is, uAs=use_As, uIhl=use_Ihl, uAd=use_Ad, uAp=use_Ap, ucapA=use_capA, urandIs=use_randIs, urandAs=use_randAs, urandN=randN)
                } else if (genotype == "wt") {
                    trials[[i]]<-runAndCalcWT_mod2(uIs=use_Is, uAs=use_As, uIhl=use_Ihl, uAd=use_Ad, uAp=use_Ap, ucapA=use_capA, urandIs=use_randIs, urandAs=use_randAs, urandN=randN)
                } else
                {
                    unknownGenoErr()
                }
                
            }
            else {
                #print("skipping!")
            }
            i.Ihl<-i.Ihl+1
            i<-i+1
        }
        i.Is<-i.Is+1
    }

    trials <- trials[lengths(trials)!=0]

    # for rbindlist
    library('data.table')

    #turn our list of data.frames into one big data.frame
    trials_df <- rbindlist(trials)

    fixnames <- function(filename) {
        if (!file.exists(csvfilename)) {
            return(NA)
            # col.names=NA add a blank name to fix alignment
            # from: https://stackoverflow.com/questions/17890777/why-are-exported-row-names-always-off-by-1-in-r
        } else {
            return(FALSE) 
        }
    }

    write.table(trials_df, csvfilename, 
                sep = ",", col.names = fixnames(csvfilename), append = T)


    source('plotParamScan.R')
    if (genotype == "tdc-1") {
        plotParamScan(uName=OUTPUT_FILENAME, 
                      #Is=20, As=2, Ihl=10, Ad=1.75, Ap=0.20 n=6
                      fitIhl=10.5, fitIs=19, fitAp=0.20, fitAd=1.75, fitAs=2, 
                      geno=genotype)
    } else if (genotype == "flp-11 nlp-7") {
        plotParamScan(uName=OUTPUT_FILENAME, 
                      fitIhl=3.5, fitIs=6.5, fitAp=0.125, fitAd=0.1667, fitAs=2,
                      geno=genotype)
    } else if (genotype == "wt") {
        plotParamScan(uName=OUTPUT_FILENAME, 
                      fitIhl=4, fitIs=7, fitAp=0.125, fitAd=1.75, fitAs=2,
                      geno=genotype)
    } else {
        unknownGenoErr()
    }
}

#================================================
#=============WT SCAN ACTIVATOR==================
#================================================

# scans values of As and Ad in WT model

gridscanner_AsAd("gridscan_wt_model2_As_Ad", As_seq=seq(0.25,6,.25), Ad=seq(0.25,6,0.25), Is=7, Ihl=4, Ap=0.125, genotype="wt", randN=6)

#================================================
#=============WT SCAN INHIBITOR==================
#================================================

# scans values of Is and Ihl in WT model

#gridscanner_IsIhl("gridscan_wt_model2_Is_Ihl", As=2, Ad=1.75, Is_seq=seq(0,14,0.5), Ihl_seq=seq(0,14,0.5), Ap=0.125, genotype="wt")

#================================================
#=============TDC-1 SCAN INHIBITOR===============
#================================================

# scans values of Is and Ihl in tdc-1 model

#gridscanner_IsIhl("gridscan_tdc1_model2_Is_Ihl", As=2, Ad=1.75, Is_seq=seq(0,33,1), Ihl_seq=seq(0,17,0.5), Ap=0.20, genotype="tdc-1")

#================================================
#==========FLP-11 NLP-7 SCAN INHIBITOR===========
#================================================

# scans values of Is and Ihl in flp-11 nlp-7 model

#gridscanner_IsIhl("gridscan_flp11nlp7_model2_Is_Ihl", As=2, Ad=1.75, Is_seq=seq(0,14,0.5), Ihl_seq=seq(0,17,0.5), Ap=0.125, genotype="flp-11 nlp-7")

#================================================
#============NLP-3 SCAN ACTIVATOR================
#================================================

# scans values of As and Ad in nlp-3 model

#gridscanner_AsAd("gridscan_nlp3_model2_As_Ad", As_seq=seq(0.1,3,0.1), Ad=c(seq(0.25,8.25,0.25)), Is=7, Ihl=4, Ap=0.06, genotype="nlp-3", randN=1)

#================================================
#============TPH-1 SCAN ACTIVATOR================
#================================================

# scans values of As and Ad in tph-1 model

#gridscanner_AsAd("gridscan_tph1_model2_As_Ad", As_seq=seq(0.1,3,.1), Ad=seq(0.25,8,0.25), Is=7, Ihl=4, Ap=0.0525, genotype="tph-1", randN=6)

