mfuzz.plotpdf <- function (eset, cl, mfrow = c(1, 1), colo, min.mem = 0, time.labels, 
    filename="fmuzzout.pdf",xlab="Time",ylab="Expression changes") 
{
    clusterindex <- cl[[3]]
    memship <- cl[[4]]
    memship[memship < min.mem] <- -1
    colorindex <- integer(dim(exprs(eset))[[1]])
    if (missing(colo)) {
        colo <- c("#FF8F00", "#FFA700", "#FFBF00", "#FFD700", 
            "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00", "#AFFF00", 
            "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00", 
            "#20FF00", "#08FF00", "#00FF10", "#00FF28", "#00FF40", 
            "#00FF58", "#00FF70", "#00FF87", "#00FF9F", "#00FFB7", 
            "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", "#00CFFF", 
            "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF", 
            "#0040FF", "#0028FF", "#0010FF", "#0800FF", "#2000FF", 
            "#3800FF", "#5000FF", "#6800FF", "#8000FF", "#9700FF", 
            "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF", 
            "#FF00D7", "#FF00BF", "#FF00A7", "#FF008F", "#FF0078", 
            "#FF0060", "#FF0048", "#FF0030", "#FF0018")
    }
    colorseq <- seq(0, 1, length = length(colo))
    for (j in 1:max(clusterindex)) {
        tmp <- exprs(eset)[clusterindex == j, ]
        tmpmem <- memship[clusterindex == j, j]
        if (((j - 1)%%(mfrow[1] * mfrow[2])) == 0) {
            pdf(filename, height=3*mfrow[1],width=3*mfrow[2])
            par(mfrow = mfrow, cex=0.5)
            if (sum(clusterindex == j) == 0) {
                ymin <- -1
                ymax <- +1
            }
            else {
                ymin <- min(tmp)
                ymax <- max(tmp)
            }
            plot.default(x = NA, xlim = c(1, dim(exprs(eset))[[2]]), 
                ylim = c(ymin, ymax), xlab = xlab, ylab = ylab, 
                main = paste("Cluster", j), axes = FALSE)
            if (missing(time.labels)) {
                axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]))
                axis(2)
            }
            else {
                axis(1, 1:dim(exprs(eset))[[2]], time.labels)
                axis(2)
            }
        }
        else {
            if (sum(clusterindex == j) == 0) {
                ymin <- -1
                ymax <- +1
            }
            else {
                ymin <- min(tmp)
                ymax <- max(tmp)
            }
            plot.default(x = NA, xlim = c(1, dim(exprs(eset))[[2]]), 
                ylim = c(ymin, ymax), xlab = xlab, ylab = ylab, 
                main = paste("Cluster", j), axes = FALSE)
            if (missing(time.labels)) {
                axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]))
                axis(2)
            }
            else {
                axis(1, 1:dim(exprs(eset))[[2]], time.labels)
                axis(2)
            }
        }
        if (!(sum(clusterindex == j) == 0)) {
            for (jj in 1:(length(colorseq) - 1)) {
                tmpcol <- (tmpmem >= colorseq[jj] & tmpmem <= 
                  colorseq[jj + 1])
                if (sum(tmpcol) > 0) {
                  tmpind <- which(tmpcol)
                  for (k in 1:length(tmpind)) {
                    lines(tmp[tmpind[k], ], col = colo[jj])
                  }
                }
            }
        }
    }
dev.off()
}


