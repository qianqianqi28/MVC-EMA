c_err=function(x){            #convexity error see Eq.36 in Weltje(1997)
    
    prop<-length(x[x<0])/length(x)
    mean_dist<-mean(x[x<0])^2
    
    con_err<-log10(prop)+log10(mean_dist)
    list(con_err=con_err)
}