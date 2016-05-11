#classifier
NBLDA <- function (x, y, xte = NULL, disperhat, beta = 1,  type = c("mle", 
    "deseq", "quantile"), prior = NULL) 
{

    if (is.null(prior)) 
        prior <- rep(1/length(unique(y)), length(unique(y)))

    null.out <- NullModel(x, type = type)
    ns <- null.out$n
    nste <- NullModelTest(null.out, x, xte, type = type)$nste

    uniq <- sort(unique(y))

    #ds <- GetDnb(ns, x, y, beta)
    ds <- matrix(1, nrow=length(uniq), ncol=ncol(x))
    for(k in 1:length(uniq)){
    	a <- colSums(x[y==uniq[k],])+beta
      	b <- colSums(ns[y==uniq[k],])+beta
      	ds[k,] <- a/b
    }

    phihat <- as.numeric(disperhat)
    discriminant <- matrix(NA, nrow = nrow(xte), ncol = length(uniq))

    for (k in 1:length(uniq)) {
        for(l in 1:nrow(xte))   {
             dstar = ds[k,]
             part2=1+nste[l,]*dstar*phihat 
             part1=dstar/part2 
             discriminant[l, k]<- sum(xte[l,]*log(part1))-sum((1/phihat)*log(part2))+log(prior[k])
         }
    }
    save <- list(ns = ns, nste = nste, ds = ds, discriminant = discriminant, 
	ytehat = uniq[apply(discriminant, 1, which.max)], x = x, y = y, xte = xte, type = type)
    return(save)
}
#generate the simulation data
SimulateDataSet <- function (n, p, K, param, sdsignal, drate) 
{
    if (n < 4 * K) 
        stop("We require n to be at least 4*K.")
    q0 <- rexp(p, rate = 1/25)
    
    isDE <-  runif(p)< drate    #runif(p)
    classk <- matrix(NA, nrow = K, ncol = p)
    for (k in 1:K) {
        lfc <- rnorm(p, sd = sdsignal)
        classk[k, ] <- ifelse(isDE, q0 * exp(lfc), q0)
    }
    truesf <- runif(n) * 2 + 0.2
    truesfte <- runif(n) * 2 + 0.2
    conds <- sample(c(rep(1:K, 4), sample(1:K, n - 4 * K, replace = TRUE)))
    condste <- sample(c(rep(1:K, 4), sample(1:K, n - 4 * K, replace = TRUE)))
    x <- xte <- matrix(NA, nrow = n, ncol = p)
    for (i in 1:n) {
        for (k in 1:K) {
            if (conds[i] == k) 
                x[i, ] <- rnbinom(p, mu = truesf[i] * classk[k, 
                  ], size = param)
            if (condste[i] == k) 
                xte[i, ] <- rnbinom(p, mu = truesfte[i] * classk[k, 
                  ], size = param)
        }
    }
    rm <- apply(x, 2, sum) == 0
    return(list(x = x[, !rm], xte = xte[, !rm], y = conds, yte = condste, 
        truesf = truesf, truesfte = truesfte, dkg=classk))
}

NBLDA_Toy_Example <- function () 
{
	library("PoiClaClu")
	library("sSeq")

	#------ Generating simulated data -------------

	truephi=20

	dat <- SimulateDataSet(n=8,p=20,sdsignal=5,K=2,param=1/truephi,drate=0.8)

	#------ Computing the dispersion using sSeq -------------

	X=t(dat$x)
	tt=getT(X,sizeFactors=rep(1,ncol(X)))$target  

	rM = rowMeans(X);
	rV = rowVars(X);

	disp = (rV - rM)/rM^2
	disp0=numeric()
	for(i in 1:length(disp)){
		 a1=0;a2=disp[i];
		 disp0[i]=max(a1,a2)
	}
	disperhat = getAdjustDisp(disp0,shrinkTarget=tt)$adj

	#------ Classification using NBLDA -------------

	out=NBLDA(dat$x,dat$y,dat$xte, disperhat)

	out$ytehat
}
