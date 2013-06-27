grandistrib <-
  function(x,scale="fine",xlab="Stations",ylab="Percentage")
  {
    
    if (scale=="fine")
      
      Descript <- granstat(x,aggr=F)$sedim[-c(1:5),]
    if (scale=="large")
      Descript <- granstat(x,aggr=F)$sedim[c(2:5),]
    
    Descript=Descript[nrow(Descript):1,]
    layout(matrix(c(1,1,1,1,1,2), nrow = 1), widths = c(0.7, 0.3))
    par(mar = c(5, 4, 4, 2) + 0.1)
    barplot(as.matrix(Descript),col=rainbow(nrow(Descript)),font.lab=2,xlab=xlab,ylab=ylab,cex.lab=1.3)
    par(mar = c(1, 1, 1, 1) + 0.1)
    plot(1,1,type="n",axes=F,xlab="",ylab="")
    legend("center",legend=rev(row.names(Descript)),fill=rev(rainbow(nrow(Descript))))
    
    
  }
