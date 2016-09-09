#https://arundurvasula.wordpress.com/2016/02/24/haplotype-networks-in-r/

c <- ape::read.nexus.data("data", format='fasta')
d <- as.DNAbin(c)
e <- dist.dna(d)
h <- pegas::haplotype(d)
h <- sort(h, what = "label")
(net <- pegas::haploNet(h))
ind.hap<-with(
  stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, pop=rownames(d)[values])
)
plot(net, size=attr(net, "freq"), scale.ratio=0.2, pie=ind.hap)
legend(-8, 0, colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=19, ncol=2)


msat<-read.genepop("twophase_2.txt",ncode=3L)
pairwise.fst(msat)
str(msat)
library(hierfstat)
pairwise.fst(msat)
fstat(msat)
Fst(msat)
library("pegas")
Fst(msat)

Fst(as.loci(msat))
pca1 <- dudi.pca(msat,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
pca1
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
s.label(pca1$li)
add.scatter.eig(pca1$eig[1:20], 3,1,2)
s.class(pca1$li, pop(msat))
add.scatter.eig(pca1$eig[1:20], 3,1,2)
s.class(pca1$li, pop(msats),xax=1,yax=3, col=transp(col,.6), axesell=FALSE
        cstar=0, cpoint=3, grid=FALSE)

latdists<-read.table("distances.txt")
latdists<-dist(latdists)

msats_strata<-genind2gtypes(msats)
structest<-structureRun(msats_strata, k.range=c(1:5),label="test",noadmix=F,pop.prior="locprior")
structurePlot(structest$test.k3.r1$q.mat)

bleh<-genindtogenpop(msat)
bleh2<-dist.genpop(bleh,method=2)

latdists<-read.table("distances.txt")
latdists<-dist(latdists)
