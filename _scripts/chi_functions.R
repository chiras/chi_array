########################################################################################
######## LIMMA ANALYSIS ################################################################
########################################################################################
textlevel = leveltext("Chi-Limma function","up",textlevel)
chi_limma <- function(eset,design,contrast.matrix,name){
textlevel = leveltext(paste("Limma Analysis (",name,")",sep=""),"up",textlevel)

fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
fit3 <- fit2

textlevel = leveltext("Mapping Annotation Data","up",textlevel)
fit2$genes$Symbol=getSYMBOL(fit2$genes$ID, "hgu133plus2")
fit2$genes$GeneName <- unlist(mget(fit$genes$ID, hgu133plus2GENENAME))
fit2$genes$Path <- paste(mget(fit$genes$ID, hgu133plus2PATH))
fit2$genes$EG <- getEG(fit2$genes$ID, "hgu133plus2")

textlevel = leveltext("Calculating Statistics","keep",textlevel)
#F.stat <- FStat(fit3)
#p.value <- pf(F.stat,df1=attr(F.stat,"df1"),df2=attr(F.stat,"df2"),lower.tail=FALSE)
topnumber = length(featureNames(eset))
tab <- topTable(fit2, adjust = "fdr", number=topnumber) # FDR corrected p-values nach dem cut off
#head(tab)


textlevel = leveltext("Writing results to disk","keep",textlevel)
write.table(tab, file=paste(PATHoutp,"/",analysis,'/results_',name,'_top',topnumber,'-',normalisierung,'.txt', sep=""),sep=";")
return(tab)
textlevel = leveltext("","down",textlevel)

}

########################################################################################
######## FILTERING DATA ################################################################
########################################################################################


textlevel = leveltext("Chi-Filtering function","keep",textlevel)
chi_filter <- function(eset,tab,IQR_filter,logFC_filter,AFFX_filter,IQR_intens_above,IQR_intens_probes,IQR_greater_than,logFC_threshold){
textlevel = leveltext("Filtering data","up",textlevel)


myfilter = ""
filtered = tab$ID
number = length(tab$ID)

textlevel = leveltext("","up",textlevel)

if (logFC_filter){
	number_orig=number


	logFC_f1 <- tab[tab$logFC > logFC_threshold | tab$logFC < -logFC_threshold,1]
	logFC_f2 <- logFC_f1  %in% featureNames(eset)
	filterOne = logFC_f1[logFC_f2]
	filtered = intersect(filterOne,filtered)
	myfilter = paste(myfilter,"_logFC",sep ="") 
	number = length(filtered)
	textlevel = leveltext(paste("logFC Filter: ", number_orig, "->",number," probesets",sep=""),"keep",textlevel)

}


if (AFFX_filter){
	number_orig=number

	affx_probes=grep("AFFX",tab$ID)

	filterTwo <- tab$ID[tab$adj.P.Val<min(tab$adj.P.Val[affx_probes])]
	filtered = intersect(filterTwo,filtered)

	myfilter = paste(myfilter,"_AFFX",sep ="") 
	number = length(filtered)
	textlevel = leveltext(paste("AFFX Filter: ", number_orig, "->",number," probesets",sep=""),"keep",textlevel)

}

### OPTION 2:  Filtern nach IQR log ratio
if (IQR_filter){
	number_orig=number

	IQR_f1 			<- pOverA(IQR_intens_probes, log2(IQR_intens_above)) 				
	IQR_f2 			<- function(x) (IQR(x) > IQR_greater_than) 		
	ff <- filterfun(IQR_f1,IQR_f2) ### Filterfunction
	filterThree <- genefilter(eset, ff) 
	filterThree  = names(filterThree)[filterThree]
	filtered = intersect(filterThree,filtered)

	myfilter = paste(myfilter,"_IQR",sep ="") 
	number = length(filtered)
	textlevel = leveltext(paste("IQR Filter: ", number_orig, "->",number," probesets",sep=""),"keep",textlevel)

}


eset <- eset[filtered, ] 
return(eset)
}

########################################################################################
######## Venn Preparation ##############################################################
########################################################################################

textlevel = leveltext("Chi-Venn function","keep",textlevel)
chi_eset2sigs <- function (eset, tab.orig){
	
	textlevel = leveltext("Defining significant rows in tabs","up",textlevel)
	inters= !tab.orig$ID %in% featureNames(eset)
	#tab.tmp = 2*(tab.orig$t>0)-1 NOT VERIFIED
	#tab.tmp[inters] = 0
	return(tab.orig[!inters,])
}

chi_eset2bins <- function (eset, tab.orig){
	textlevel = leveltext("Returning significant hits","up",textlevel)

	inters= !tab.orig$ID %in% featureNames(eset)
	tab.tmp = tab.orig$t
	tab.tmp[tab.orig$t>0] = 1 
	tab.tmp[tab.orig$t<0] = -1 
	tab.tmp[tab.orig$t==0] = 0 

	tab.tmp[inters] = 0
	return(tab.tmp)
}


chi_venn <- function(setlistPLUS,setlistMINUS){
if (length(setlistPLUS)>1){

OLlistPLUS <- overLapper(setlist=setlistPLUS, sep="_", type="vennsets")
OLlistMINUS <- overLapper(setlist=setlistMINUS, sep="_", type="vennsets")

pdf(file=paste(PATHoutp,"/",analysis,'/MFA-Plots/Venn_',paste(names(setlistPLUS), collapse = "-"),'.pdf', sep=""))
	counts <- list(sapply(OLlistMINUS$Venn_List, length), sapply(OLlistPLUS$Venn_List, length));
	vennPlot(counts=counts, mysub="Top: Up-regulated; Bottom: Down-regulated")
dev.off()

}}

chi_overlap <- function(setlistPLUS,setlistMINUS){
	overlap = list()
	overlap$both_up = list()
	overlap$both_down = list()
	overlap$contrary = list()

	for (i in 1:length(setlistPLUS)){
	for (j in i:length(setlistPLUS)){
		if(i != j){
			namecomp=paste(names(setlistPLUS)[i],"--",names(setlistPLUS)[j],sep="")
			overlap$both_up[[namecomp]]=setlistPLUS[[i]][setlistPLUS[[i]] %in% setlistPLUS[[j]]]
			if(length(overlap$both_up[[namecomp]])==0){overlap$both_up[[namecomp]]=NA}
		}
	}}
	
	for (i in 1:length(setlistMINUS)){
	for (j in i:length(setlistMINUS)){
		if(i != j){
			namecomp=paste(names(setlistMINUS)[i],"--",names(setlistMINUS)[j],sep="")
			overlap$both_down[[namecomp]]=setlistMINUS[[i]][setlistMINUS[[i]] %in% setlistMINUS[[j]]]
			if(length(overlap$both_down[[namecomp]])==0){overlap$both_down[[namecomp]]=NA}
		}
	}}

	for (i in 1:length(setlistPLUS)){
	for (j in i:length(setlistMINUS)){
		if(i != j){
			namecomp=paste(names(setlistPLUS)[i],"--",names(setlistMINUS)[j],sep="")
			overlap$contrary[[namecomp]]=setlistPLUS[[i]][setlistPLUS[[i]] %in% setlistMINUS[[j]]]
			if(length(overlap$contrary[[namecomp]])==0){overlap$contrary[[namecomp]]=NA}
		}
	}}
	return(overlap)
}


########################################################################################
######## Heat Preparation #(OLD)########################################################
########################################################################################

	
textlevel = leveltext("Chi-Heatmap preparation function","keep",textlevel)
Probesets2Heatmap <- function(eset.norm, VectorOfProbesets){ 
eset2plot <- eset.norm[VectorOfProbesets, ] # Expressionsset für diese Probes generieren

i=1
while (file.exists(paste("R-Output/Paths/path_ind",i,".pdf",sep=""))){
	i=i+1
}
pdf(file=paste("R-Output/Paths/path_ind",i,".pdf",sep=""),  height=length(featureNames(eset2plot))/2+2,width=8)
	
heatmap.2(exprs(eset2plot), labRow=getSYMBOL(featureNames(eset2plot), "hgu133plus2"), Rowv = F, Colv = F,dendrogram="none",lmat=rbind( c(0, 4,0,0), c(2,1,1,0 ), c(0,3,0,0) ), lwid=c(0.5, 2, 4,0.5 ), lhei= c(2,length(featureNames(eset2plot))/2,0.1), density.info="none", scale="row", key = T, symkey = T, trace ="none")
	 # # hier wird die eigentliche Heatmap erstellt!
mtext("           Individual Heatmap Selection", side=2, line=-2, adj=0.0, cex=1, col="black", outer=TRUE)
dev.off()

}

textlevel = leveltext("Chi-Pathway to Heatmap function","keep",textlevel)
PathSigHeat <- function(eset_sig, pathwaynr){
	
cat(noquote(paste("Heatmap: KEGG", pathwaynr,": Obtaining pathway associated probesets: ", sep="")))
eset.sig.norm.path <- eset.sig.norm[grep(pathwaynr,kegg)]

symbs<- mget(featureNames(eset.sig.norm.path), hgu133plus2SYMBOL, ifnotfound = NA)

symbs<- unlist(symbs)

used_symbs=c()
cat(noquote(paste(length(symbs),"/",length(featureNames(eset.sig.norm.path)), " of ",length(kegg),"\n", sep="")))

cat(noquote(paste("Heatmap: KEGG", pathwaynr,": Getting Gene Symbols\n", sep="")))
for (i in 1:length(symbs)){
	xxnum=0;
	while (sum(used_symbs == paste(symbs[i], xxnum, sep="-"))){
		xxnum= xxnum +1;
	}
	symbs[i]=paste(symbs[i], xxnum, sep="-")	
	used_symbs = c(used_symbs, symbs[i])
}

head(symbs)

cat(noquote(paste("Heatmap: KEGG", pathwaynr,": Defining new Matrix\n", sep="")))
mat <- exprs(eset.sig.norm.path)
cat(noquote(paste("Heatmap: KEGG", pathwaynr,": Renaming Matrix\n",sep="")))
rownames(mat) <- ifelse(!is.na(symbs), as.vector(symbs), featureNames(eset.sig.norm.path))

cat(noquote(paste("Heatmap: KEGG", pathwaynr,": Plotting to PDF\n", sep="")))
pdf(file=paste('R-Output/Paths/path_top',topnumber,'-KEGG', pathwaynr,'-',auswertung,'-',normalisierung,'-',myfilter,'.pdf', sep=""), height=length(rownames(mat))/2+2,width=8)

heatmap.2(mat[order(rownames(mat)), ] , Rowv = F, Colv = F,dendrogram="none", lmat=rbind( c(0, 4,0,0), c(2,1,1,0 ), c(0,3,0,0) ), lwid=c(0.5, 2, 4,0.5 ), lhei= c(2,length(rownames(mat))/2,0.1),density.info="none", scale="row", key = T, symkey = T, trace ="none")
mtext(paste("           ",as.character(mget(pathwaynr, KEGGPATHID2NAME))), side=2, line=-2, adj=0.0, cex=1, col="black", outer=TRUE)

dev.off()

}

########################################################################################
######## Heatmap ##############################################################
########################################################################################


textlevel = leveltext("Chi-Heatmap function","keep",textlevel)
heatmap.chi <-function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
distfun = dist, hclustfun = hclust, dendrogram = c("both",
"row", "column", "none"), symm = FALSE, scale = c("none",
"row", "column"), na.rm = TRUE, revC = identical(Colv,
"Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) ||
scale != "none", col = "heat.colors", colsep, rowsep,
sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1,
notecol = "cyan", na.color = par("bg"), trace = c("column",
"row", "both", "none"), tracecol = "cyan", hline = median(breaks),
vline = median(breaks), linecol = tracecol, margins = c(5,
5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr),
cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL,
key = TRUE, keysize = 1.5, density.info = c("histogram",
"density", "none"), denscol = tracecol, symkey = min(x <
0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL,
xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL,
...)
{
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
    "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
    "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
        c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
            dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
            dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
        c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
            dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
            dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
        stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm)
        x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
        stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
    (1:nr)[rowInd]
    else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
    (1:nc)[colInd]
    else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) <
    1) {
        if (missing(col) || is.function(col))
        breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks)
        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
        length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
    col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || length(ColSideColors) !=
            nc)
            stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] +
            1)
            lhei <- c(lhei[1], 0.2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || length(RowSideColors) !=
            nr)
            stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) -
            1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], 0.2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
        ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
    c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col,
    breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
    retval$rowDendrogram <- ddr
    if (exists("ddc"))
    retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
        col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
    cex.axis = cexCol)
    if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
    cex.axis = cexRow)
    if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
    eval(substitute(add.expr))
    if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0,
    length(csep)), xright = csep + 0.5 + sepwidth[1],
    ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1,
    col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) +
    1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) +
    1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1,
    col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
    col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    #  else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    #  else plot.new()
    if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
        xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row")
        mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column")
        mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, "Value", line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    #   else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
    high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}

########################################################################################
######## Interaktion-plot ##############################################################
########################################################################################


textlevel = leveltext("Chi-Interaction Plot function","keep",textlevel)
InteractKeggPlot <- function(top_table, pathName){ 

	top= top_table
	pathName = pathName
	x <- hgu133plus2ENTREZID
	top$ENTREZ <- unlist(as.list(x[top$ID]))
	top <- top[!is.na(top$ENTREZ),] 
	top <- top[!duplicated(top$ENTREZ),] 
	tg1 <- top #[top$adj.P.Val < 0.05,]
	DE_tg1 <- tg1$logFC
	names(DE_tg1) <- as.vector(tg1$ENTREZ)	
	ALL_tg1 <- top$ENTREZ

	tmp <- paste("../KEGGdata/hsa", pathName,".xml",sep="")
	
	if (!file.exists(tmp)){download.file(paste("http://www.genome.jp/kegg-bin/download?entry=hsa", pathName,"&format=kgml",sep=""), tmp)}
	
	pdf(file=paste('R-Output/Paths/path_inter',topnumber,'-KEGG', pathName,'-',auswertung,'-',normalisierung,'-',myfilter,'.pdf', sep=""), width=20, height=20)
	g <- parseKGML2Graph(tmp)
	
	deKID <- translateGeneID2KEGGID(names(DE_tg1))
	allKID <- translateGeneID2KEGGID(ALL_tg1)
	
	isDiffExp <- nodes(g) %in% deKID
	# sprintf("%2.2f%% genes differentially-expressed", mean(isDiffExp)*100)
	
	g.old = g
	g <- subGraph(nodes(g)[isDiffExp], g)
	
	ar <- 20
	cols <- heatPalette(ar) # colorRampPalette(brewer.pal(6, "RdBu"))(ar)
	colsO <- heatPalette(2) # colorRampPalette(brewer.pal(6, "RdBu"))(ar)
	
	logfcs <- DE_tg1[match(nodes(g), deKID)]
	names(logfcs) <- nodes(g)
	logfcs[is.na(logfcs)] <- 0
#	incol <- round((logfcs+2)*5);
	incol <- round((logfcs+min(logfcs)*-1)*5);

	outcol <-logfcs
	outcol[outcol>0]<-2
	outcol[outcol<0]<-1
	incol[incol>ar] <- ar
	undetected <- !nodes(g) %in% allKID
		
	logcol <- cols[incol]; 
	logcol[logfcs==0] <- "darkgrey"; 
	logcol[undetected] <- "yellow";
	
	logcolO <- colsO[outcol]; 
	
	names(logcol) <- names(logfcs)
	nA <- makeNodeAttrs(g, fillcolor=logcol, color= logcolO, label="", width=10, height=1.2,shape="triangle")
	par(mar=c(3,5,0,5), mgp=c(0,0,0))
	
	gGeneID <- translateKEGGID2GeneID(nodes(g))
	gSymbol <-  sapply(gGeneID, function(x) mget(x, org.Hs.egSYMBOL, ifnotfound=NA)[[1]])
		
	#gnA$width <- makeAttr(g,list("0.8"=toprbccName))
	
	layout(mat=matrix(c(rep(1,8),2), ncol=1, byrow=TRUE))
	gnA <- makeNodeAttrs(g, fillcolor=logcolO, lwd=3, label=gSymbol, fontsize=1,fixedsize=T)
	plot(g, "twopi", nodeAttrs=gnA)
	dev.off()
}

########################################################################################
######## GO stats ##############################################################
########################################################################################


textlevel = leveltext("Chi-GOstats function","keep",textlevel)
chi_gostats <- function(analysisname, eset.sig.norm,eset.norm){

textlevel = leveltext("Obtaining Gene Clusters and Universe","up",textlevel)

geneCluster=as.vector(na.omit(getEG(featureNames(eset.sig.norm),"hgu133plus2")))
geneCluster2=as.vector(na.omit(getSYMBOL(featureNames(eset.sig.norm),"hgu133plus2")))
geneCluster3=featureNames(eset.sig.norm)

geneUniverse=as.vector(na.omit(getEG(featureNames(eset.norm),"hgu133plus2"))) 
geneUniverse2=as.vector(na.omit(getSYMBOL(featureNames(eset.norm),"hgu133plus2"))) 
geneUniverse3=featureNames(eset.norm)

gos = c("KEGG","MF","BP","CC")
gosigs=list()

topnumber=length(featureNames(eset.sig.norm))

outfile = paste(PATHoutp,"/",analysis,'/GO_',analysisname,'_','top',topnumber,'-PATH-',normalisierung,'.txt', sep="")

textlevel = leveltext("Creating output file","keep",textlevel)

sink(outfile)
writeLines(paste("","GOid","Type","Name","Pvalue","OddsRatio","ExpCount","Count","Size","GOGenesNumber","SigGenesNumber","SigGenesLinCorr","GOGenesSymbol","GOGenesEG","SigGenesSymbol","SigGenesEG",sep=";"), con=paste(PATHoutp,"/",analysis,'/GO_',analysisname,'_','top',topnumber,'-PATH-',normalisierung,'.txt', sep=""))	
sink()


for (id in 1:4){
	textlevel = leveltext(paste("Getting overrepresented groups for: ",gos[id], sep=""),"keep",textlevel)
	textlevel = leveltext("Calculation","up",textlevel)

	if (gos[id] != "KEGG"){

		paramsGOunder <- new("GOHyperGParams", geneIds = geneCluster, universeGeneIds = geneUniverse, annotation = "hgu133plus2", ontology = gos[id], pvalueCutoff = 1, conditional = FALSE, testDirection = "under")
		paramsGOover <- new("GOHyperGParams", geneIds = geneCluster, universeGeneIds = geneUniverse, annotation = "hgu133plus2", ontology = gos[id], pvalueCutoff = 1, conditional = FALSE, testDirection = "over")
		tryCatch(hgOverGO <- hyperGTest(paramsGOover),error = function(e) {print('error GO')})
		#tryCatch(hgUnderGO <- hyperGTest(paramsGOunder),error = function(e) {print('error GO')})
		#htmlReport(hgUnderGO,file=paste('GOunder_',gos[id],'_top',topnumber,'-PATH-',auswertung,'-',normalisierung,'.html', sep=""))
		#htmlReport(hgOverGO,file=paste('GOover_',gos[id],'_top',topnumber,'-PATH-',auswertung,'-',normalisierung,'.html', sep=""))

	}else{

		paramsGOunder <- new("KEGGHyperGParams", geneIds = geneCluster, universeGeneIds = geneUniverse, annotation = "hgu133plus2",pvalueCutoff = 1, testDirection = "under")
		paramsGOover <- new("KEGGHyperGParams", geneIds = geneCluster, universeGeneIds = geneUniverse, annotation = "hgu133plus2",pvalueCutoff = 1, testDirection = "over")
		tryCatch(hgOverGO <- hyperGTest(paramsGOover),error = function(e) {print('error KEGG')})
		#tryCatch(hgUnderKEGG <- hyperGTest(paramsKEGGunder),error = function(e) {print('error KEGG')})
		#htmlReport(hgOverKEGG,file=paste('KEGGover_top',topnumber,'-PATH-',auswertung,'-',normalisierung,'-',myfilter,'.html', sep=""))
		#htmlReport(hgUnderKEGG,file=paste('KEGGunder_top',topnumber,'-PATH-',auswertung,'-',normalisierung,'-',myfilter,'.html', sep=""))
	}
	
	textlevel = leveltext("Summarizing","keep",textlevel)
	gosigs[[gos[id]]]=summary(hgOverGO)[summary(hgOverGO)$Pvalue<0.05,][[1]]
	output = gos[id]

	listofGOs = gosigs[[output]]

if (length(listofGOs)>0){
	textlevel = leveltext("Parsing","keep",textlevel)
	textlevel = leveltext("","up",textlevel)

for (i in 1:length(listofGOs)){ 
	
analyzeGO=F;
if (gos[id] == "KEGG" & length(unlist(mget(listofGOs[i],hgu133plus2PATH2PROBE, ifnotfound=NA))[is.na(unlist(mget(listofGOs[i],hgu133plus2PATH2PROBE, ifnotfound=NA)))==F])>0){analyzeGO<-T}
if (gos[id] != "KEGG" & length(unlist(mget(listofGOs[i],hgu133plus2GO2PROBE  , ifnotfound=NA))[is.na(unlist(mget(listofGOs[i],hgu133plus2GO2PROBE  , ifnotfound=NA)))==F])>0){analyzeGO<-T}

if(analyzeGO){
genefunction=list();

	textlevel = leveltext(paste("GO-Set: ",gos[id] ," - ", i, " (",length(listofGOs),") ", sep=""),"keep",textlevel)

#cat(noquote(paste("GOstats: Working on GO-Set: ",gos[id] ," - ", i, " (",length(listofGOs),") ", sep="")))
#cat(noquote("[1] Collecting Data"))
myinfo = summary(hgOverGO)[summary(hgOverGO)[,1]==listofGOs[i],]

#cat(noquote(" [2] Parsing Data"))
genefunction$GOid					= myinfo[,1]	
genefunction$Type					= names(myinfo)[1]	
genefunction$Name 					= myinfo$Term
genefunction$Pvalue 				= myinfo$Pvalue
genefunction$OddsRatio 				= myinfo$OddsRatio
genefunction$ExpCount 				= myinfo$ExpCount
genefunction$Count 					= myinfo$Count
genefunction$Size 					= myinfo$Size

if(gos[id] != "KEGG"){
	genefunction$GOGenesNumber 			= length(unique(as.vector(unlist(getSYMBOL(unlist(mget(listofGOs[i],hgu133plus2GO2PROBE), use.names=F),"hgu133plus2"), use.names=F))))
	tmpGOGenesSymbol 					= unique(as.vector(unlist(getSYMBOL(unlist(mget(listofGOs[i],hgu133plus2GO2PROBE), use.names=F),"hgu133plus2"), use.names=F)))
	tmpGOGenesEG 						= unique(as.vector(unlist(getEG(unlist(mget(listofGOs[i],hgu133plus2GO2PROBE), use.names=F),"hgu133plus2"), use.names=F)))
	#genefunction$AllGenesChipNumber 	= length(unique(intersect(tmpGOGenesSymbol, geneUniverse2)))
	genefunction$SigGenesNumber 		= length(unique(intersect(tmpGOGenesSymbol, geneCluster2)))
	genefunction$SigGenesCorrected 		= (genefunction$Size/genefunction$GOGenesNumber)*genefunction$SigGenesNumber
	genefunction$GOGenesSymbol 			= paste(unique(as.vector(unlist(getSYMBOL(unlist(mget(listofGOs[i],hgu133plus2GO2PROBE), use.names=F),"hgu133plus2"), use.names=F))), collapse=", ")
	genefunction$GOGenesEG 				= paste(unique(as.vector(unlist(getEG(unlist(mget(listofGOs[i],hgu133plus2GO2PROBE), use.names=F),"hgu133plus2"), use.names=F))), collapse=", ")
	#genefunction$AllGenesChipSymbol 	= paste(unique(intersect(tmpGOGenesSymbol, geneUniverse2)), collapse=", ")
	#genefunction$AllGenesChipEG 		= paste(unique(intersect(tmpGOGenesEG, geneUniverse)), collapse=", ")
	genefunction$SigGenesSymbol 		= paste(unique(intersect(tmpGOGenesSymbol, geneCluster2)), collapse=", ")
	genefunction$SigGenesEG				= paste(unique(intersect(tmpGOGenesEG, geneCluster)), collapse=", ")
}else{
	genefunction$GOGenesNumber 			= length(unique(as.vector(unlist(getSYMBOL(unlist(mget(listofGOs[i], hgu133plus2PATH2PROBE), use.names=F),"hgu133plus2"), use.names=F))))
	tmpGOGenesSymbol 					= unique(as.vector(unlist(getSYMBOL(unlist(mget(listofGOs[i], hgu133plus2PATH2PROBE), use.names=F),"hgu133plus2"), use.names=F)))
	tmpGOGenesEG 						= unique(as.vector(unlist(getEG(unlist(mget(listofGOs[i], hgu133plus2PATH2PROBE), use.names=F),"hgu133plus2"), use.names=F)))
	#genefunction$AllGenesChipNumber 	= length(unique(intersect(tmpGOGenesSymbol, geneUniverse2)))
	genefunction$SigGenesNumber 		= length(unique(intersect(tmpGOGenesSymbol, geneCluster2)))
	genefunction$SigGenesCorrected 		= (genefunction$Size/genefunction$GOGenesNumber)*genefunction$SigGenesNumber
	genefunction$GOGenesSymbol 			= paste(unique(as.vector(unlist(getSYMBOL(unlist(mget(listofGOs[i], hgu133plus2PATH2PROBE), use.names=F),"hgu133plus2"), use.names=F))), collapse=", ")
	genefunction$GOGenesEG 				= paste(unique(as.vector(unlist(getEG(unlist(mget(listofGOs[i], hgu133plus2PATH2PROBE), use.names=F),"hgu133plus2"), use.names=F))), collapse=", ")
	#genefunction$AllGenesChipSymbol 	= paste(unique(intersect(tmpGOGenesSymbol, geneUniverse2)), collapse=", ")
	#genefunction$AllGenesChipEG 		= paste(unique(intersect(tmpGOGenesEG, geneUniverse)), collapse=", ")
	genefunction$SigGenesSymbol 		= paste(unique(intersect(tmpGOGenesSymbol, geneCluster2)), collapse=", ")
	genefunction$SigGenesEG				= paste(unique(intersect(tmpGOGenesEG, geneCluster)), collapse=", ")	
}

#cat(noquote(" [3] Writing Data\n"))
#cat(noquote(paste(c("",t(as.data.frame(genefunction))),collapse=";")),file=outfile, append=T, sep="\n")

}else{
	textlevel = leveltext(paste("GO-Set excluded: ",gos[id] ," - ", i, " (",length(listofGOs),")", sep=""),"keep",textlevel)
}# end if (analyzeGO)
if (devel_mode && i > 2) {break}

}# end for length(listofGOs))
textlevel = leveltext("","down",textlevel)

}# end if (length(listofGOs)>0)

textlevel = leveltext("Making Barplot","keep",textlevel)

cat(noquote(paste(c("",t(as.data.frame(genefunction))),collapse=";")),file=outfile, append=T, sep="\n")

sig.hgOverGO=summary(hgOverGO)[summary(hgOverGO)$Pvalue<0.1,]

sig.hgOverGO =sig.hgOverGO[order(sig.hgOverGO$Pvalue, decreasing=T),]

pdf(file=paste(PATHoutp,"/",analysis,'/SA-Plots/GO_',analysisname,'_',gos[id],'_top',topnumber,'-',normalisierung,'.pdf', sep=""), height=4+length(sig.hgOverGO$Count)/10,width=10)
 	par(mar=c(3,20,0.1,0.2))
	barplot(rbind(sig.hgOverGO$ExpCount,sig.hgOverGO$Count), col=c("gray80","gray20"), beside=T, horiz=T, names.arg=paste(sig.hgOverGO$KEGGID,sig.hgOverGO$Term,"- p =", round(sig.hgOverGO$Pvalue, digits=4)), las=2, cex.names=0.5)


 	for (i in 1:(max(sig.hgOverGO$Count)/5)){
 		abline(v=i*5, lty=3)
 	}
 	par(xpd=NA)
 	legend(-15,0,legend=c("Count", "Expected"), fill=c("gray20","gray80"))
dev.off()

textlevel = leveltext("","down",textlevel)
}# end for MF/BP/CC/KEGG
}# end function analyzeGOstats


########################################################################################
######## PATHWAY ANALYSIS ##############################################################
########################################################################################

textlevel = leveltext("Chi-Genes2Pathway function\n","keep",textlevel)
chi_genes2Pathways<- function(name, tab){
 	tab2 = tab[1,]
 	tab2[1,] = c(rep("NA",length(tab2[1,])))
	rownames(tab2) = "Start"

textlevel = leveltext("Parsing","up",textlevel)
for (i in 1:(length(tab$Path))){
 #	print(i);
 	if (is.na(tab$Path[i]) | tab$Path[i]  == "NA"){
 		tab2 =rbind(tab2,tab[i,])
 	}else{
 		pathxways= gregexpr("[0-9]+", tab$Path[i])
 		for (j in 1:length(attr(pathxways[[1]],"match.length"))){
 			bla=substr(tab$Path[i],pathxways[[1]][j],pathxways[[1]][j]+attr(pathxways[[1]],"match.length")[j]-1)
 			tab2 =rbind(tab2,tab[i,])
 			tab2$Path[length(tab2$Path)]=paste(bla,"-",unlist(mget(bla, KEGGPATHID2NAME)))
 			}
 			
 	}#endelse
 	}#endif

topnumber=length(tab[,1])

textlevel = leveltext("Writing results to disk","keep",textlevel)
write.table(tab2, file=paste(PATHoutp,"/",analysis,'/pwt_',name,'_top',topnumber,'-',normalisierung,'.txt', sep=""),sep=";")

textlevel = leveltext("Plotting abundance diagram","keep",textlevel)
tab3=table(tab2$Path[tab2$Path!="NA"])[order(table(tab2$Path[tab2$Path!="NA"]),decreasing = F)]
pdf(file=paste(PATHoutp,"/",analysis,'/SA-Plots/PW_',name,'_top',topnumber,'-',normalisierung,'.pdf', sep=""), height=2+length(tab3[tab3>2])/10,width=10)
 	par(mar=c(3,20,0.1,0.2))
 	barplot(tab3[tab3>2], horiz=T, las=2, cex.names=0.75)
 	for (i in 1:(max(tab3)/5)){
 		abline(v=i*5, lty=3)
 	}
dev.off()

textlevel = leveltext("Plotting abundance diagram with only pathways > 10","keep",textlevel)
tab3=table(tab2$Path[tab2$Path!="NA"])[order(table(tab2$Path[tab2$Path!="NA"]),decreasing = F)]
pdf(file=paste(PATHoutp,"/",analysis,'/SA-Plots/PW_',name,'_top',topnumber,'-',normalisierung,'_10plus.pdf', sep=""), height=2+length(tab3[tab3>9])/10,width=10)
 	par(mar=c(3,20,0.1,0.2))
 	barplot(tab3[tab3>9], horiz=T, las=2, cex.names=0.75)
 	for (i in 1:(max(tab3)/5)){
 		abline(v=i*5, lty=3)
 	}
dev.off()

}



