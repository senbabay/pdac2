library(shiny)
###################################################################
vars = readRDS("PDAC_Affy_82_assayData_gene_variances.rds") # variance before standardization
zmat1 = readRDS("dataparts/PDAC_Affy_82_assayData_scaled_part1.rds") # standardized data matrix
zmat2 = readRDS("dataparts/PDAC_Affy_82_assayData_scaled_part2.rds")
zmat3 = readRDS("dataparts/PDAC_Affy_82_assayData_scaled_part3.rds")
zmat = rbind(zmat1,zmat2,zmat3)
pheno = read.csv("PDAC_Affy_82_phenoData.txt",sep="\t",header=T,stringsAsFactors = F)
features = readRDS("PDAC_Affy_82_featureData_small.rds")
gsymbols = features$Gene.Symbol
# Short-term vs long-term stats
ls.pval = readRDS("PDAC_Affy_82_long_vs_short_pval.rds")
ls.boolean = readRDS("PDAC_Affy_82_long_vs_short_boolean.rds")
### Unsupervised clustering
# Threshold vector for variances
lx = seq(3,50,1) / 10
pve = len = rep(NA,length(lx))

# Number of genes passing cutoff, and percent variance explained at that step
for(i in 1:length(lx)){
  #print(i)
  SDCUT = lx[i]
  len[i] = length(which(vars > SDCUT))
  
  zmat2 = zmat[which(vars > SDCUT),]
  pc = prcomp(t(zmat2))
  x3 = pc$x[,1:3]
  pve[i] = sum(apply(x3,2,sd)^2) / sum(apply(pc$x,2,sd)^2)
}#end for i

respCol = rep("magenta",nrow(pheno))
respCol[which(pheno$Cohort == "long")] = "turquoise4"
mylab = pheno$Well

###################################################################
shinyServer(function(input, output) { 
  
  varcut <- reactive({
    input$animation
  })
  
  output$vardensity <- renderPlot({
    vc = varcut()
    dx = density(vars)
    
    par(mar=c(4,4,3,1)+0.1)
    plot(dx,log="x",lwd=3,col=4,main="Distribution of probe set variances",xlab="Probe set variance across 30 samples",cex.lab=1.5,cex.axis=1.5)
    abline(v=vc,col=2,lwd=3)
    mtext(text="Probe set var cutoff",side=3,line=0,at=vc,cex=1,col=2)
    
    ty = floor(max(dx$y))-0.5
    w = which(lx==vc)
    text(vc,ty,paste(len[w],"\nprobe sets",sep=""),pos=4,offset=0.3,col=2,font=2)
  })
  
  output$percentVar <- renderPlot({
    vc = varcut()
    w = which(lx==vc)
    par(mar=c(4,4,3,1)+0.1)
    plot(lx[1:w],100*pve[1:w],pch=16,type="p",xlab="Probe set variance cutoff",ylab="% variance explained",cex.lab=1.5,col=4,
         main="Percent variance explained by first 3 PCs",xlim=range(lx),ylim=c(100*min(pve),100*max(pve)),cex.axis=1.5,cex=1.5)
  })
  
  getScoreMat <- reactive({
    vc = varcut()
    zmat2 = zmat[which(vars > vc),]
    pc = prcomp(t(zmat2))
    pc$x
  })
  
  output$pc1pc2 <- renderPlot({
    x = getScoreMat()
    par(mar=c(4,4,3,1)+0.1)
    plot(x[,1],x[,2],pch=16,xlab="PC1",ylab="PC2",cex.lab=1.5,col=respCol,
         main="PC1 vs PC2",cex.axis=1.5,cex=1.5)
    text(x[,1],x[,2],labels=mylab,cex=0.6,pos=1,offset=0.5)
    legend("topleft",horiz=T,pch=16,legend=c("LT","ST"),col=c("turquoise4","magenta"),pt.cex=1.5,cex=0.8)
  })
  
  output$pc1pc3 <- renderPlot({
    x = getScoreMat()
    par(mar=c(4,4,3,1)+0.1)
    plot(x[,1],x[,3],pch=16,xlab="PC1",ylab="PC3",cex.lab=1.5,col=respCol,
         main="PC1 vs PC3",cex.axis=1.5,cex=1.5)
    text(x[,1],x[,3],labels=mylab,cex=0.6,pos=1,offset=0.5)
    legend("topleft",horiz=T,pch=16,legend=c("LT","ST"),col=c("turquoise4","magenta"),pt.cex=1.5,cex=0.8)
  })
  
  getTable <- reactive({
    vc = varcut()
    w = which(vars > vc)
    
    # significance between LT and ST
    pw = ls.pval[w]; hiLT = ls.boolean[w]
    
    # variance
    tmp = vars[w]
    # data frame
    df = data.frame(PROBE.SET=names(tmp),GENE.SYMBOL=gsymbols[w],VARIANCE=as.numeric(tmp),HIGH.in.LT=hiLT,PVAL=signif(pw,3))
    # order according to variance
    ordx = order(df$VARIANCE,decreasing=T)
    df = df[ordx,]
    #uniqGene = sapply(strsplit(as.character(df[,1]),"//"),"[",c(T,F))
    #url = paste0("http://www.genecards.org/cgi-bin/carddisp.pl?gene=",uniqGene)
    #df[,1] = paste0("<a href='",url,"' target='_blank'>",df[,1],"</a>")
    df
  })
  
  output$geneTable <- DT::renderDataTable({
    DT::datatable(getTable(),
                  escape=FALSE,
                  rownames = FALSE,
                  filter="top",
                  options = list(pageLength = 20),
                  caption = 'Probe sets sorted from highest to lowest variance')
  })
  
  output$slideHeatmap <- renderUI({
    vc = varcut()
    sliderInput("slideHeatmap", "Probe set variance cutoff:", 0.3, 5, vc,step = 0.1)
  })  
  
  output$heatmap <- renderImage({
    vc = input$slideHeatmap
    filename <- normalizePath(file.path('../shinyHeatmaps',paste0('cutoff',vc, '.png')))
    #filename <- normalizePath(file.path(paste0('cutoff',vc, '.png')))
    # Return a list containing the filename and alt text
    list(src = filename,
         width = 900,
         alt = paste("Cutoff", vc))
  }, deleteFile = FALSE)
  
  
  output$legend <- renderImage({
    filename <- normalizePath(file.path('./',"legend_small.png"))
    list(src = filename,
         width = 300,
         height=100,
         alt = "Heatmap legend")
  }, deleteFile = FALSE)
  
  
  output$numGene <- renderText({
    w = which(lx==input$slideHeatmap)
    paste0(len[w],"\nprobe sets")
  })
}
)




