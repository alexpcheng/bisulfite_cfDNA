# Measuring mixing parameters
mp_function<-function(reference_and_samples, other, nt, sum_to_one){
  num.tissues<-nt
  A<-as.data.frame(reference_and_samples[,4:(num.tissues+3)])
  b_list<-as.list(reference_and_samples[,(num.tissues+4):ncol(reference_and_samples)])
  sample.name<-names(b_list)
  percentages<-lapply(b_list, mixing_parameters_function, A, other=other, sample.name, sum_to_one)
  return(percentages)
}

mixing_parameters_function<-function(b, A, other, sample_name, sum_to_one){
  #Accounting for missing tissues
  if(other){
    A$other<-1
    G<-diag(ncol(A))
    g<-rep(-1,ncol(A))
    g[ncol(A)]<-0
    G<-rbind(G, g)
    h<-rep(0,ncol(A))
    h[ncol(A)+1]<-(-1)
    E<-NULL
    f<-NULL
    sol<-lsei(A=A, B=b, G=G, H=h, E=E, F=f, type=2, fulloutput = T) #calls solve.QP from package quadprog
    sol.df<-data.frame(sol$X)
    sol.df$tissue<-rownames(sol.df)
    #sol.df<-sol.df[!rownames(sol.df) %in% "other",]
    if (sum_to_one){
      normalizing_factor<-sum(sol.df$sol.X)
      sol.df$sol.X<-sol.df$sol.X/normalizing_factor
    }
    rel_error<-(sol$solutionNorm/nrow(A))
    rel_error<-data.frame(rel_error, "rel_error")
    colnames(rel_error)<-c(sample_name, "tissue")
  }
  #Not Accounting for missing tissues
  if(!other){
    G<-diag(ncol(A))
    h<-rep(0,ncol(A))
    E<-(rep(1,ncol(A)))
    f<-1
    sol<-lsei(A=A, B=b, G=G, H=h, E=E, F=f, type=1, fulloutput = T)
    sol.df<-as.data.frame(sol$X)
    sol.df$tissue<-rownames(sol.df)
    if (sum_to_one){
      normalizing_factor<-sum(sol.df$sol.X)
      sol.df$sol.X<-sol.df$sol.X/normalizing_factor
    }
    rel_error<-(sol$solutionNorm/nrow(A))
    rel_error<-data.frame(rel_error, "rel_error")
    colnames(rel_error)<-c(sample_name, "tissue")
  }
  colnames(sol.df)[1]<-sample_name
  sol.df<-rbind(sol.df, rel_error)
  return(sol.df)
}
