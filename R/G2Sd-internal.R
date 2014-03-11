.fowa.stat <-
function(x,phi,um){
    
    folk.ward=data.frame(matrix(ncol=0,nrow=9))
    
    
    for (b in 1:dim(x)[2])
    {y=x[b]
     sum.sieve=sum(y)
     class.weight=(y*100)/sum.sieve
     cum.sum=cumsum(class.weight)[,1]
     
     if (min(cum.sum)>5) 
     {
       fowa=data.frame(rep(0,9))
       row.names(fowa)=c("Sediment","Mean.fw.um","Sd.fw.um","Skewness.fw.um","Kurtosis.fw.um","Mean.fw.phi","Sd.fw.phi","Skewness.fw.phi","Kurtosis.fw.phi")
       names(fowa)=names(x)[b]
     }
     
     if (min(cum.sum)<5)
     {
       mat.D=.percentile(x[,b],phi,um)
       
       mean.phi=(mat.D[6,1]+mat.D[4,1]+mat.D[2,1])/3
       mean.mm=(exp(log(mat.D[6,2])+log(mat.D[4,2])+log(mat.D[2,2])/3))/1000
       
       sd.phi=-(((mat.D[2,1]-mat.D[6,1])/4)+((mat.D[1,1]-mat.D[7,1])/6.6))
       sd.mm=exp(((log(mat.D[2,2])-log(mat.D[6,2]))/4)+((log(mat.D[1,2])-log(mat.D[7,2]))/6.6))
       
       skewness.phi=-(((mat.D[6,1]+mat.D[2,1]-(2*mat.D[4,1]))/(2*(mat.D[2,1]-mat.D[6,1])))+ ((mat.D[7,1]+mat.D[1,1]-(2*mat.D[4,1]))/(2*(mat.D[1,1]-mat.D[7,1]))))
       skewness.mm=-skewness.phi
       
       kurtosis.phi=(mat.D[1,1]-mat.D[7,1])/(2.44*(mat.D[3,1]-mat.D[5,1]))
       kurtosis.mm=kurtosis.phi
       
       if (mean.phi<=-5) mean.descript="Very Coarse Gravel"
       if (mean.phi>-5 & mean.phi<=-4) mean.descript="Coarse Gravel"
       if (mean.phi>-4 & mean.phi<=-3) mean.descript="Medium Gravel"
       if (mean.phi>-3 & mean.phi<=-2) mean.descript="Fine Gravel"
       if (mean.phi>-2 & mean.phi<=-1) mean.descript="Very Fine Gravel"
       if (mean.phi>-1 & mean.phi<=0) mean.descript="Very Coarse Sand"
       if (mean.phi>0 & mean.phi<=1) mean.descript="Coarse Sand"
       if (mean.phi>1 & mean.phi<=2) mean.descript="Medium Sand"
       if (mean.phi>2 & mean.phi<=3) mean.descript="Fine Sand"
       if (mean.phi>3 & mean.phi<=4) mean.descript="Very Fine Sand"
       if (mean.phi>4 & mean.phi<=5) mean.descript="Very Coarse Silt"
       if (mean.phi>5 & mean.phi<=6) mean.descript="Coarse Silt"
       if (mean.phi>6 & mean.phi<=7) mean.descript="Medium Silt"
       if (mean.phi>7 & mean.phi<=8) mean.descript="Fine Silt"
       if (mean.phi>8 & mean.phi<=9) mean.descript="Very Fine Silt"
       if (mean.phi>8) mean.descript="Clay"
       
       if (sd.phi<0.35) sorting="Very Well Sorted"
       if (sd.phi>=0.35 & sd.phi<0.5) sorting="Well Sorted"
       if (sd.phi>=0.5 & sd.phi<0.7) sorting="Moderately Well Sorted"
       if (sd.phi>=0.7 & sd.phi<1) sorting="Moderately Sorted"
       if (sd.phi>=1 & sd.phi<2) sorting="Poorly Sorted"
       if (sd.phi>=2 & sd.phi<4) sorting="Very Poorly Sorted"
       if (sd.phi>=4) sorting="Extremely Poorly Sorted"
       
       
       if (skewness.phi>=0.3) skewness.descript="Very Fine Skewed"
       if (skewness.phi<0.3 & skewness.phi>=0.1) skewness.descript="Fine Skewed"
       if (skewness.phi<0.1 & skewness.phi>-0.1) skewness.descript="Symmetrical"
       if (skewness.phi<=-0.1 & skewness.phi>-0.3) skewness.descript="Coarse Skewed"
       if (skewness.phi<=-0.3) skewness.descript="Very Coarse Skewed"
       
       if (kurtosis.phi<0.67) kurtosis.descript="Very Platykurtic"
       if (kurtosis.phi>=0.67 & kurtosis.phi<0.9) kurtosis.descript="Platykurtic"
       if (kurtosis.phi>=0.9 & kurtosis.phi<=1.11) kurtosis.descript="Mesokurtic"
       if (kurtosis.phi>1.11 & kurtosis.phi<=1.5) kurtosis.descript="Leptokurtic"
       if (kurtosis.phi>1.5 & kurtosis.phi<=3) kurtosis.descript="Very Leptokurtic"
       if (kurtosis.phi>3) kurtosis.descript="Extremely Leptokurtic"
       
       .sedim.descript=paste(mean.descript,sorting,skewness.descript,kurtosis.descript,sep=",")
       
       result.fw.phi=data.frame(c(round(mean.phi,3),round(sd.phi,3),round(skewness.phi,3),round(kurtosis.phi,3)))
       names(result.fw.phi)=names(x)[b]
       
       result.fw.mm=data.frame(c(round(mean.mm,3),round(sd.mm,3),round(skewness.mm,3),round(kurtosis.mm,3)))
       names(result.fw.mm)=names(x)[b]
       
       fowa=data.frame(rbind(.sedim.descript,result.fw.mm,result.fw.phi))
       row.names(fowa)=c("Sediment","Mean.fw.um","Sd.fw.um","Skewness.fw.um","Kurtosis.fw.um","Mean.fw.phi","Sd.fw.phi","Skewness.fw.phi","Kurtosis.fw.phi")
       names(fowa)=names(x)[b]
     }
     folk.ward=cbind(folk.ward,fowa)
    }
    folk.ward
  }
.G2Sd_web <-
function(){
  runApp(appDir=paste0(.libPaths(),"/G2Sd/extdata"))
}
.grancompat <-
function(x)
{
  x <- as.data.frame(x)
  n.sieve <- nrow(x)
  n.sample <- ncol(x)
  if (!is.data.frame(x)) 
    stop("dataframe expected.")
  if (any(x < 0))
    stop("negative entries in dataframe.")
  if (any(x > 300))
    warning("Some high values are present.", call. = FALSE,immediate.=TRUE)
  if (n.sieve!=29)
  {
    cat("Compatibility progress.... \n \n")
    
    ref_sieve=c(25000,20000,16000,12500,10000,8000,6300,5000,4000,2500,
                2000,1600,1250,1000,800,630,500,400,315,250,200,160,125,
                100,80,63,50,40,0)      
    
    init_df <- as.data.frame(matrix(data=0,ncol=n.sample,nrow=length(ref_sieve)));colnames(init_df) <- colnames(x)
    row.names(init_df) <-ref_sieve
    
   
    if (any(is.na(pmatch(row.names(x),ref_sieve))))
      stop("Incorrect sieve values.")
    
    else 
   
    {for (sieve in row.names(x))
      init_df[sieve,] <- x[sieve,]}
    
  }
  else init_df <- x
  return(init_df)
}
.index.sedim <-
function(x,phi,um){
    x=as.data.frame(x)
    INDEX=data.frame(matrix(ncol=0,nrow=9))
    for (b in 1:dim(x)[2])
    {
      mat.D=.percentile(x[,b],phi,um)
      index=data.frame(matrix(ncol=1,nrow=9))
      row.names(index)=c("D10(um)","D50(um)","D90(um)","D90/D10","D90-D10","D75/D25","D75-D25","Trask(So)","Krumbein(Qd)")
      names(index)=names(x)[b]
      index[1,1]=round(mat.D[9,2],3)
      index[2,1]=round(mat.D[4,2],3)
      index[3,1]=round(mat.D[8,2],3)
      index[4,1]=round(mat.D[8,2]/mat.D[9,2],3)
      index[5,1]=round(mat.D[8,2]-mat.D[9,2],3)
      index[6,1]=round(mat.D[3,2]/mat.D[5,2],3)
      index[7,1]=round(mat.D[3,2]-mat.D[5,2],3)
      index[8,1]=round(sqrt(mat.D[3,2]/mat.D[5,2]),3)
      index[9,1]=round((mat.D[5,1]-mat.D[3,1])/2,3)
      
      
      INDEX=cbind(INDEX,index)
    }
    return(INDEX)
  }
.mgrep <-
function(mpattern,x,FUN="grep")
{
  select=NULL
  for (i in 1:length(mpattern))
  {
    if (FUN=="grep")
      values=grep(mpattern[i],x)
    if (FUN=="which")
      values=which(x==mpattern[i])  
    select=c(select,values)
  }
  return(sort(select))
}
.mode.sedim <-
function(x,um){
    
    x=as.data.frame(x)
    sum.sieve=apply(x,2,sum)
    
    MODE=data.frame(matrix(ncol=0,nrow=5))
    for (b in 1:dim(x)[2])
    {
      
      
      class.weight=(x[,b]*100)/sum.sieve[b]            
      tab.mod=cbind(um,class.weight)
      if (pmatch(0,um)!=0) tab.mod=tab.mod[-pmatch(0,um),]
      
      plot(tab.mod[,1],tab.mod[,2],type="b",lwd=3,xlab="Particule size (microns)",ylab="Pourcentage (%)",xaxt="n",log="x")
      a=identify(tab.mod,plot=FALSE,n=4)
      
      mod=data.frame(tab.mod[a,1])
      names(mod)=names(x)[b]
      row.names(mod)=tab.mod[a,1]
      
      if (dim(mod)[1]==1) mod.descript="1 Mode" else mod.descript=paste(dim(mod)[1],"Modes")
      
      
      MODE.sedim=data.frame(matrix(ncol=1,nrow=4))
      
      for (i in 1:dim(mod)[1])
        MODE.sedim[i,]=mod[i,]
      MODE.sedim=rbind(mod.descript,MODE.sedim)
      names(MODE.sedim)=names(x)[b]
      row.names(MODE.sedim)[1]="Nb Mode"
      MODE.sedim
      MODE=cbind(MODE,MODE.sedim)
    }
    return(MODE)
  }
.moment.arith <-
function(x,um){
    
    
    x=as.data.frame(x)
    sum.sieve=apply(x,2,sum)
    
    arith=data.frame(matrix(ncol=0,nrow=4))
    for (b in 1:dim(x)[2])
    {
      class.weight=(x[b]*100)/sum.sieve[b]
      
      
      
      mid.point=rep(0,(length(um)))
      
      for(i in 2:length(um))
      {
        
        mid.point[i]=(um[i]+um[i-1])/2
        
      }
      
      fm=class.weight*mid.point
      mean.arith=apply(fm,2,sum)/100
      
      fmM2=class.weight*(mid.point-mean.arith)^2
      sd.arith=sqrt(apply(fmM2,2,sum)/100)
      
      fmM3=class.weight*(mid.point-mean.arith)^3
      skewness.arith=apply(fmM3,2,sum)/(100*sd.arith^3)
      
      fmM4=class.weight*(mid.point-mean.arith)^4
      kurtosis.arith=apply(fmM4,2,sum)/(100*sd.arith^4)
      
      
      moment.arit=data.frame(rbind(round(mean.arith,3),round(sd.arith,3),round(skewness.arith,3),round(kurtosis.arith,3)))
      colnames(moment.arit)=colnames(x)[b]
      arith=cbind(arith,moment.arit)
      rownames(arith)=c("mean.arith.um","sd.arith.um","skewness.arith.um","kurtosis.arith.um")
    }
    return(arith)
  }
.moment.geom <-
function(x,phi){
    
    
    x=as.data.frame(x)
    sum.sieve=apply(x,2,sum)
    
    geom=data.frame(matrix(ncol=0,nrow=4))
    for (b in 1:dim(x)[2])
    {
      class.weight=(x[b]*100)/sum.sieve[b]
      
      mid.point=rep(0,(length(phi)))
      
      for(i in 2:length(phi))
      {
        
        mid.point[i]=(phi[i]+phi[i-1])/2
        
      }
      
      
      logm=log10(2^(-mid.point)*1000)
      flogm=class.weight*logm
      mean.geom=10^(apply(flogm,2,sum)/100)
      
      fmM2=class.weight*(logm-log10(mean.geom))^2
      sd.geom=10^(sqrt(apply(fmM2,2,sum)/100))
      
      fmM3=class.weight*(logm-log10(mean.geom))^3
      skewness.geom=(apply(fmM3,2,sum)/(100*log10(sd.geom)^3))
      
      fmM4=class.weight*(logm-log10(mean.geom))^4
      kurtosis.geom=(apply(fmM4,2,sum)/(100*log10(sd.geom)^4))
      kurtosis3.geom=kurtosis.geom
      
      moment.geo=as.data.frame(rbind(round(mean.geom,3),round(sd.geom,3),round(skewness.geom,3),round(kurtosis.geom,3)))
      names(moment.geo)=names(x)[b]
      geom=cbind(geom,moment.geo)
      rownames(geom)=c("mean.geom.um","sd.geom.um","skewness.geom.um","kurtosis.geom.um")
    }
    return(geom)
  }
.northarrow <-
function(loc,size,bearing=0,cols,letter_dist=1,cex=1,...) {
    # checking arguments
    if(missing(loc)) stop("loc is missing")
    if(missing(size)) stop("size is missing")
    # default colors are white and black
    if(missing(cols)) cols <- rep(c("white","black"),8)
    # calculating coordinates of polygons
    radii <- rep(size/c(1,4,2,4),4)
    x <- radii[(0:15)+1]*cos((0:15)*pi/8+bearing)+loc[1]
    y <- radii[(0:15)+1]*sin((0:15)*pi/8+bearing)+loc[2]
    # drawing polygons
    for (i in 1:15) {
      x1 <- c(x[i],x[i+1],loc[1])
      y1 <- c(y[i],y[i+1],loc[2])
      polygon(x1,y1,col=cols[i])
    }
    # drawing the last polygon
    polygon(c(x[16],x[1],loc[1]),c(y[16],y[1],loc[2]),col=cols[16])
    # drawing letters
    b <- c("E","N","W","S")
    for (i in 0:3) text((size+letter_dist*par("cxy")[1])*cos(bearing+i*pi/2)+loc[1],
                        (size+letter_dist*par("cxy")[2])*sin(bearing+i*pi/2)+loc[2],b[i+1],
                        cex=cex)
  }
.percentile <-
function(x,phi,um){
    
    x=as.numeric(x)
    sum.sieve=sum(x)
    class.weight=(x*100)/sum.sieve
    cum.sum=as.numeric(cumsum(class.weight))
    D=c(5,16,25,50,75,84,95,10,90)
    
    if (min(cum.sum)>5)
    {warning("D5 can't be calculated : 
             Folk-Ward statistics can't be computed.", call. = FALSE,immediate.=TRUE)
     mat.D=data.frame(matrix(ncol=2,nrow=9,0))
     row.names(mat.D)=100-D
     names(mat.D)=c("Phi","um")}
    
    if (min(cum.sum)<5){
      
      class.weight.PHI=cbind(cum.sum,phi,um)         
      mat.D=data.frame(matrix(ncol=2,nrow=9))
      row.names(mat.D)=100-D
      names(mat.D)=c("Phi","um")
      
      for (i in 1:9)
      {
        greaterpercent=subset(class.weight.PHI,class.weight.PHI[,1]>D[i])
        greaterphi=subset(greaterpercent,greaterpercent[,1]==min(greaterpercent[,1]),select=-1)
        greaterphi=as.numeric(subset(greaterphi,greaterphi[,2]==max(greaterphi[,2]),select=1))
        greaterpercent=min(greaterpercent[,1])
        lesspercent=subset(class.weight.PHI,class.weight.PHI[,1]<D[i])
        lessphi=subset(lesspercent,lesspercent[,1]==max(lesspercent[,1]),select=-1)
        lessphi=as.numeric(subset(lessphi,lessphi[,2]==min(lessphi[,2]),select=1))
        lesspercent=max(lesspercent[,1])
        
        ifelse  (dim(subset(class.weight.PHI,class.weight.PHI[,1]==D[i]))[1]==0,
{ratio1=(D[i]-lesspercent)/(greaterpercent-lesspercent)
 ratio2=(greaterphi-lessphi)*ratio1
 phi=lessphi+ratio2
 um=1/(2^phi)*1000},
{phi=as.numeric(subset(class.weight.PHI,class.weight.PHI[,1]==D[i],2))
 um=as.numeric(subset(class.weight.PHI,class.weight.PHI[,1]==D[i],3))*1000}) 
        
        
        result=c(phi,um)
        mat.D[i,]=result
      }
    }
    return(mat.D)
    
  }
.Random.seed <-
c(403L, 406L, -1838796285L, -1731966943L, -48438512L, -2083618578L, 
1814667417L, 588711211L, 1857371162L, 1686212316L, -12179745L, 
-1440198475L, -1329621428L, 1279949474L, 312943101L, -2114606345L, 
1306999838L, 622491352L, 1348577851L, -920008743L, 2064158056L, 
-939524970L, -838141423L, -248941629L, 2125696978L, -690206108L, 
596412839L, 1449797373L, -419519756L, 1499014618L, 1856370981L, 
-617291713L, -184688698L, 2036614224L, -655529869L, -1218401135L, 
-483570496L, 1670839966L, 1775948393L, -914468197L, 373694378L, 
1966975948L, 1914889871L, 54965957L, -980910468L, 1806947026L, 
-1159604659L, 1959688967L, -426066194L, -81855544L, -575173429L, 
2146727209L, -343337608L, 2104986278L, -1054643007L, 289901907L, 
-150599166L, 237146932L, 782237879L, 579709101L, -2010982972L, 
-819386902L, 687618645L, -185482257L, 1085170678L, -25963232L, 
1998846691L, 1817959297L, 248659440L, 544105614L, 268250105L, 
-1308112501L, -1768589894L, 341499068L, 561655807L, 323669333L, 
-1487311252L, 870984386L, -1362416355L, -34888809L, -282315010L, 
-89332552L, -1876142437L, 1784870137L, 127239752L, 168518006L, 
-887404367L, 237341923L, -1867531342L, -876452860L, 1522181959L, 
549027357L, -1945234284L, 1724445754L, -1599123835L, 1310830431L, 
-1155665946L, 1557185136L, -1740429933L, -440028623L, -352735008L, 
702278014L, 572754953L, 91053115L, 1966203402L, -1437734420L, 
898608751L, -1134100059L, -17622692L, -345832142L, -1889518547L, 
93508327L, -1497573362L, 1005963752L, 1639844843L, 304690313L, 
-1514000360L, -451973306L, -132305503L, -1306565581L, 1132538274L, 
-393048940L, -1178096873L, -1096686323L, -1144990684L, -2074770422L, 
-1918035595L, -620150577L, -443190314L, 1192061824L, 1562763203L, 
-860956319L, -881741872L, -1200052434L, 577682137L, -2119329557L, 
1536846170L, 436450076L, 580052511L, -1403011467L, -1400857844L, 
1156814562L, 81558461L, 848600119L, 692276958L, -2085344744L, 
-1718454917L, 1128836377L, -1413910104L, 1641514326L, -1045989935L, 
-192344701L, 347296786L, -148103260L, 1599879015L, -604281667L, 
96097332L, 1451988890L, -1623342107L, -280739969L, 811196166L, 
-1612506992L, -156389069L, 1506780881L, 739783808L, 417609054L, 
520915625L, -1463052965L, 451184490L, -1414803956L, 1121797711L, 
-822539515L, 1751321532L, 197697170L, 1216118157L, -213602873L, 
-702378322L, -2010028280L, -1401326325L, 1511916521L, -212454216L, 
1882694502L, 479315969L, 1198823059L, 1141975362L, 139986676L, 
-90369289L, -2028764435L, 87727492L, -1331293270L, -1629476587L, 
-1084922961L, 831929270L, 1677435360L, 912048419L, 745535041L, 
1551218992L, 1546619726L, -1037102663L, 1982532555L, -2109243270L, 
-1904287364L, -619220545L, -1351102827L, -974031700L, -1338237310L, 
-240863651L, 1776965975L, -1329118146L, -383661960L, -1337879717L, 
-518648135L, 375735048L, -1025246794L, 1536984817L, -322039517L, 
-2032909454L, -1137743420L, -735506041L, -1595217059L, -649616044L, 
475121530L, -837179707L, 1768531999L, 506810278L, 1662366000L, 
785984467L, 475197937L, -514287840L, 622550736L, -732032862L, 
579785768L, -1455188724L, -405147716L, 1967665458L, -1001950784L, 
1062922804L, 1912144840L, -867241798L, 1072351632L, 1360357868L, 
-1480263084L, -268093294L, -745089248L, 1688269260L, -1423493664L, 
-383080766L, -843178968L, -1016431956L, 1589978300L, -1830243854L, 
-594624208L, -255886172L, -2021740872L, 625667418L, 1488428432L, 
181484412L, -610684652L, -725019838L, 1254308608L, 1166305212L, 
81684304L, -1896504542L, 133780808L, 70414860L, -693397412L, 
-181926702L, -964468224L, 797726036L, -2111520280L, 895151098L, 
2117144784L, 1151184172L, 122880756L, 1879859602L, 986698848L, 
-778990324L, 269562592L, 1257017698L, 508574376L, 445129292L, 
873724284L, 2038037810L, 2015055984L, 2030546756L, -17062440L, 
668957818L, 1573928496L, -110193988L, -488134828L, -645970878L, 
-1508524640L, -1280955716L, -1612721712L, 1692632418L, -1602088984L, 
857882572L, -419629316L, -302229006L, 1021055104L, 1830857972L, 
-680194424L, 1314357050L, 2065128272L, 1740861036L, 905267348L, 
-637593262L, 452859616L, -1349657908L, 1394213024L, 1990069122L, 
-1852710168L, 43729644L, 1120316476L, 1783282290L, 1086456944L, 
-1123189212L, -2114347976L, 1485768730L, 1755806416L, 390093884L, 
995828884L, 407346306L, 1984143296L, 1861120316L, 1607989904L, 
-2129267678L, 1738216456L, 2045470860L, -1807195108L, 1404550802L, 
-515226496L, -638758892L, -677798616L, -607019462L, 1426869456L, 
-417831828L, 1920217780L, 1533452306L, 1497307104L, 1236146124L, 
-481686048L, 1733461538L, 838943016L, -1761765108L, 1196385852L, 
-486575886L, -474549520L, -1822538940L, 965796440L, -1783215558L, 
131890928L, -735873860L, 1049122324L, -946227006L, 1210157152L, 
1359697276L, -150842160L, -771351518L, 1981039272L, -293230708L, 
-114970180L, -1121070030L, -2115020992L, 546531508L, -1487586488L, 
-1536386374L, 992597776L, 715342060L, -1703095212L, 955899666L, 
-649634912L, -358645812L, 1081719648L, 19419330L, 551265832L, 
738436524L, 556457788L, -661387406L, -629466320L, -442088924L, 
-1151805384L, -1195776806L, -2120226416L, 1656030972L, 316141460L, 
-144774334L, -674143744L, -2063877572L, 2047737424L, 1958902434L, 
-777069752L, 1857675788L, -1570591524L, -1698932142L, -846133120L, 
-1642259628L, -96224792L, -668579078L, -1739884720L, 773667884L, 
-1721748108L, -221280622L, 2010135136L, -565525876L, 413072224L, 
929489634L, 68829352L, 2145910988L, -258553348L, 1036423090L, 
-132360720L, 892061124L, 2052661080L, -888050822L, -1308977744L, 
-1532126788L, 717105236L, -1738888254L, -174030048L, 2056945596L, 
-2058237360L, 10780130L, -977006232L, 248548684L, 882919676L, 
375808754L, 784406272L, -1946139020L, 384886792L, -42911302L, 
497653712L, -2062951956L, 909332372L, 212890322L, 18902624L, 
-468122292L, 793957664L, -1483124094L, -1411700248L, -1960671252L, 
-1557989828L, -1411043726L, 1895113584L, -925210332L, 389853752L, 
-1773272678L, 210504144L, 1633164732L, 752298132L, 68320386L, 
-844118592L, 1037583164L, -443485040L, 2097095458L, 1776922120L, 
1117246050L, 372577527L, 285462337L, 1032285638L, -745685244L, 
-334136651L, -642697825L, -1570846216L, 90458070L, 2009420387L, 
1661766181L, -778784318L, -1316342096L, 57159385L, -1863927733L, 
-608901444L, -905804566L, 370430239L, -1347788375L, 137086782L, 
2042573356L, 1838768253L, -956677433L, 877657008L, -719639090L, 
-1708187909L, 232613213L, -574542198L, -1261343576L, 597276689L, 
2068023555L, -1013488636L, 1080822834L, 1783161991L, 1480306897L, 
-412523946L, -2113385932L, 124154181L, -2105914481L, -88772856L, 
-1877157594L, -458356269L, 722101045L, -453610862L, 2012915296L, 
-791896119L, 86494971L, 1481314508L, 27855194L, -1800510385L, 
1068099865L, -67169490L, 1300036924L, -2060653907L, 749135191L, 
-750054240L, -1280237058L, -816837429L, 1396247245L, 487286266L, 
-430373320L, 1546800097L, -2012900205L, 1955925684L, -1899793470L, 
-1211004969L, -1090992095L, -966082586L, -1004152092L, 801730069L, 
1267185215L, -692902440L, -985183242L, -1306999805L, 2054350661L, 
-70093982L, 616570064L, 1372080121L, 1732004651L, 129368220L, 
2118065866L, -2009593921L, -218717431L, -2105877474L, 133785036L, 
-939004003L, 1745779239L, 1753700432L, -48414482L, -1093268069L, 
950399997L, -228274454L, 876616072L, -1186367695L, 2104207267L, 
-980160604L, 1453378130L, 179347879L, 148624881L, -2108808906L, 
-2144847404L, -653701275L, -366739601L, 1052787112L, -1725016314L, 
1782977075L, 259870869L, 2035871986L, 1353312576L, -911829335L, 
2113154843L, -784664340L, -128822022L, -163209937L, 2062125369L, 
861354830L, 1545011612L, 827972365L, -598069897L, 1574513024L, 
2146463710L, 1555683499L, 1394510765L, 1994059802L, 1254287064L, 
1157949249L, 733304819L, -1524131308L, -143899742L, 1304771767L, 
1639224961L, 490073734L, -310323516L, 839489525L, -58704289L, 
-1170794312L, -794714730L, 522506147L, 1681933157L, -350150014L, 
-587993104L, -475039335L, -611715189L, 436126332L, 1497649322L, 
-545656097L, 466731113L, -1043514754L, 1661862508L, 898563773L, 
-216939641L, 216606576L, -1615770098L, 542462779L, -1633054691L, 
-941242166L, 49222632L, 816193745L, -269255485L, 1347847108L, 
1688250482L, -1719305913L, -1299700463L, -1463000938L, -411638284L, 
564754309L, -1574328625L, 18428872L, -1370004890L, 136176420L
)
.sedim.descript <-
function(x,um){
    
    
    um=as.numeric(um)
    
    x=as.data.frame(x)
    sum.sieve=apply(x,2,sum)
    sediment=data.frame(matrix(ncol=0,nrow=13))
    for (b in 1:dim(x)[2])
    {
      
      class.weight=(x[,b])*100/sum.sieve[b]
      class.weight.um=cbind(class.weight,um)
      
      
      seuil.sedim=c(63000,31500,2^c(4:-3)*1000,63,40,NA)
      class.sedim=c("boulder","vcgravel","cgravel","mgravel","fgravel","vfgravel","vcsand","csand","msand","fsand","vfsand","vcsilt","silt")
      sedim=data.frame(cbind(seuil.sedim,class.sedim),stringsAsFactors = FALSE)
      sedim[,1]=as.numeric(sedim[,1])
      
      
      result=data.frame(matrix(nrow=dim(sedim)[1],ncol=dim(sedim)[1]))
      result[,1]=sedim[,1]
      names(result)=c("Sedim","Pourcentage")
      sedim.result=0
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]>=sedim[1,1])
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[1,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[1,1] & class.weight.um[,2]>=(sedim[2,1]))
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[2,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[2,1] & class.weight.um[,2]>=(sedim[3,1]))
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[3,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[3,1] & class.weight.um[,2]>=(sedim[4,1]))
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[4,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[4,1] & class.weight.um[,2]>=(sedim[5,1]))
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[5,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[5,1] & class.weight.um[,2]>=(sedim[6,1]))
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[6,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[6,1] & class.weight.um[,2]>=(sedim[7,1]))
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[7,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[7,1] & class.weight.um[,2]>=(sedim[8,1]))
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[8,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[8,1] & class.weight.um[,2]>=(sedim[9,1]))
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[9,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[9,1] & class.weight.um[,2]>=(sedim[10,1]))
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[10,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[10,1] & class.weight.um[,2]>=(sedim[11,1]))
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[11,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[11,1] & class.weight.um[,2]>=(sedim[12,1]))
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[12,1]=sedim.result
      
      sedim.percent=subset(class.weight.um,class.weight.um[,2]<sedim[12,1])
      sedim.result=sum(as.numeric(sedim.percent[,1]))
      result[13,1]=sedim.result
      
      result=data.frame(result[,1])
      names(result)=names(x)[b]
      row.names(result)=class.sedim
      
      sediment=cbind(sediment,round(result,3))
    }
    return(sediment)
  }
.texture.sedim <-
function(x,um){
    
    um=as.numeric(um)
    x=as.data.frame(x)
    sum.sieve=apply(x,2,sum)
    
    Texture=data.frame(matrix(ncol=0,nrow=5))
    for (b in 1:dim(x)[2])
    {
      
      class.weight=(x[,b]*100)/sum.sieve[b]           
      class.weight.um=cbind(class.weight,um)
      class.weight.um[,1] <- as.numeric(as.character(class.weight.um[,1]))
      class.weight.um[,2] <- as.numeric(as.character(class.weight.um[,2]))
      
      seuil.texture=c(63000,2000,63,0)
      class.texture=c("Boulder","Gravel","Sand","mud")
      texture=data.frame(cbind(seuil.texture,class.texture),stringsAsFactors = FALSE)
      texture[,1]=as.numeric(texture[,1])
      
      
      
      result=data.frame(matrix(nrow=dim(texture)[1],ncol=dim(texture)[2]))
      result[,1]=texture[,2]
      names(result)=c("Texture","Pourcentage")
      texture.result=0
      
      texture.percent=subset(class.weight.um,class.weight.um[,2]>=texture[1,1])
      texture.result=sum(texture.percent[,1])
      result[1,2]=texture.result
      
      texture.percent=subset(class.weight.um,class.weight.um[,2]<texture[1,1] & class.weight.um[,2]>=(texture[2,1]))
      texture.result=sum(texture.percent[,1])
      result[2,2]=texture.result 
      
      texture.percent=subset(class.weight.um,class.weight.um[,2]<texture[2,1] & class.weight.um[,2]>=(texture[3,1]))
      texture.result=sum(texture.percent[,1])
      result[3,2]=texture.result 
      
      texture.percent=subset(class.weight.um,class.weight.um[,2]<(texture[3,1]))
      texture.result=sum(texture.percent[,1])
      result[4,2]=texture.result 
      
      mud=round(result[4,2],3)
      gravel=round(result[2,2],3)
      sand=round(result[3,2],3)
      
      
{if (mud==0 & sand==0) mudsand=0
 if (mud==0 & sand>0) mudsand=10
 if (sand==0 & mud>0) mudsand=0.01
 else mudsand=sand/mud }
      
      if (mudsand>=9){
        if (gravel>80) texture="Gravel"
        if (gravel>30 & gravel<=80) texture="Sandy Gravel"
        if (gravel>5 & gravel<=30) texture="Gravelly Sand"
        if (gravel>0 & gravel<=5) texture="Slightly Gravelly Sand"
        if (gravel==0) texture="Sand"}
      
      if (mudsand>=1 & mudsand<9){
        if (gravel>80) texture="Gravel"
        if (gravel>30 & gravel<=80) texture="Muddy Sandy Gravel"
        if (gravel>5 & gravel<=30) texture="Gravelly Muddy Sand"
        if (gravel>0 & gravel<=5) texture="Slightly Gravelly Muddy Sand"
        if (gravel==0) texture="Muddy Sand"}
      
      if (mudsand>=(1/9) & mudsand<1){
        if (gravel>80) texture="Gravel"
        if (gravel>30 & gravel<=80) texture="Muddy Gravel"
        if (gravel>5 & gravel<=30) texture=" Gravelly Mud"
        if (gravel>0 & gravel<=5) texture="Slightly Gravelly Sandy Mud"
        if (gravel==0) texture="Sandy Mud"} 
      
      if (mudsand<(1/9)){
        if (gravel>80) texture="Gravel"
        if (gravel>30 & gravel<=80) texture="Muddy Gravel"
        if (gravel>5 & gravel<=30) texture=" Gravelly Mud"
        if (gravel>0 & gravel<=5) texture="Slightly Gravelly Mud"
        if (gravel==0) texture="Mud"} 
      
      row.names(result)=result[,1]
      name.texture=row.names(result)
      result=data.frame(result[,2])
      
      texture.sedim=data.frame(rbind(texture,round(result,3)))
      row.names(texture.sedim)=c("Texture",name.texture)
      names(texture.sedim)=names(x)[b]  
      texture.sedim
      Texture=cbind(Texture,texture.sedim)
    }
    
    return(Texture)
  }
