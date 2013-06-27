granplot <-
function(x,xc=1,hist=TRUE,cum=TRUE,main="",col.cum="red",col.hist="gray") 
  {

    
    x <- x[order(as.numeric(row.names(x))),]
  
    
    um=as.numeric(row.names(x))
    if (pmatch(0,um)!=0 ) um[pmatch(0,um)]="<40"
    
    
    
sum.sieve = sum(x[, xc])
class.weight = (x[, xc] * 100)/sum.sieve
class.weight.cum = cumsum(class.weight)

if (hist == TRUE & cum == TRUE) {
  par(mar = c(5, 4, 4, 4))
  barplot(x[, xc] * ((100/max(x[, xc]))), xlab = "Particule size (microns)", 
          ylab = "Weight (g)", names.arg = um, las = 2, yaxt = "n", 
          main = main, col = col.hist,cex.names=0.75,font.lab=2)

  axis(4, at = seq(0, 100, 20), labels = seq(0, 100, 20), 
       las = 2, col = col.cum, col.axis = col.cum,cex.axis=0.8)
  axis(2, at = seq(0, 100, 20), labels = seq(0, max(x[, 
                                                      xc]), max(x[, xc])/5), las = 2,cex.axis=0.8)
  mtext("(% cum) ", 4, line = 2, col = col.cum,font=2)
  lines(class.weight.cum, col = col.cum, lwd = 1)
}
if (hist == FALSE & cum == TRUE) {
  plot(class.weight.cum, type = "l", lwd = 1, xlab = "Particule size (microns)", 
       ylab = "Percentage cum.(%)", las = 2, xaxt = "n", 
       xlim = c(0, 29), main = main, col = col.cum,cex.axis=0.7,font.lab=2)

  axis(1, at = 1:29, labels = um, las = 2,cex.axis=0.7)
}
if (hist == TRUE & cum == FALSE) {
  barplot(x[, xc], xlab = "Particule size (microns)", ylab = "Weight (g)", 
          names.arg = um, las = 2, main = main, col = col.hist,cex.names=0.75,font.lab=2)
  
}
}
