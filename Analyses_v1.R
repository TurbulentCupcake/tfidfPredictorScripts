library(DECIPHER)

l <- load("/Users/sixfoot10/Desktop/Classification/AdithyaMurali/161229/loto_balanced.RData")
l <- load("/Users/sixfoot10/Desktop/Classification/AdithyaMurali/161229/loto_unbalanced.RData")

# sanity checks
all(loto_unbal$genera==loto_bal$genera) # should be TRUE
any(loto_unbal$genera==loto_unbal$predictions) # should be FALSE
any(loto_bal$genera==loto_bal$predictions) # should be FALSE

# make a tree from the selected sequence representatives
rna <- readRNAStringSet("~/Desktop/Classification/AdithyaMurali/161229/RDP_alignment_v2.fas.gz")
names(rna) <- sapply(strsplit(names(rna), ";", fixed=TRUE),
	function (x) x[length(x)])
d <- DistanceMatrix(rna[loto_bal$rdpIndex])
c <- IdClusters(d, type="dendrogram")

# compare confidences between bal and unbal
plot(loto_bal$confidences, loto_unbal$confidences)
w <- which(loto_unbal$confidences > loto_bal$confidences + 0.4)

# look at the size of genera with greater confidence in unbal than bal
t <- table(names(rna)) # genera sizes (number of representatives)
loto_unbal$genera[!(loto_unbal$genera %in% names(t))] # a few genera names are missing quotes
t[loto_unbal$predictions[w]] # only big assigned genera (> 3) have higher confidence in unbal

# look at the fold-change of predicting a bigger genus in unbal
h_unbal <- hist(t[loto_unbal$predictions], breaks=c(seq(0.5, 10, 1), 1000))
h_bal <- hist(t[loto_bal$predictions], breaks=c(seq(0.5, 10, 1), 1000))
plot(h_unbal$density/h_bal$density, ylab="Unbal/Bal", xlab="Genus size")
abline(h=1)

# look at the number of errors in unbal versus bal
sum(loto_unbal$confidences > 0.5) # total "errors" in unbal
sum(loto_unbal$confidences > 0.5 & t[loto_unbal$predictions] > 4, na.rm=TRUE) # total "errors" to big genera in unbal
sum(loto_bal$confidences > 0.5) # total "errors" in bal
sum(loto_bal$confidences > 0.5 & t[loto_bal$predictions] > 4, na.rm=TRUE) # total "errors" to big genera in bal
# (CONCLUSION: ABOUT A QUARTER OF OC ERRORS COULD BE ELIMINATED BY FIXING THE IMBALANCE ISSUE)

# plot distance versus confidence in OC errors
w <- which(loto_unbal$predictions %in% rownames(d) & loto_unbal$genera %in% rownames(d)) # only genera without quoting issues in their names
plot(d[cbind(loto_unbal$predictions[w], loto_unbal$genera[w])],
	loto_unbal$confidences[w],
	xlab="Distance to prediction",
	ylab="Confidence in prediction",
	main="Unbalanced")
w <- which(loto_bal$predictions %in% rownames(d) & loto_bal$genera %in% rownames(d)) # only genera without quoting issues in their names
plot(d[cbind(loto_bal$predictions[w], loto_bal$genera[w])],
	loto_bal$confidences[w],
	xlab="Distance to prediction",
	ylab="Confidence in prediction",
	main="Balanced") # should we investigate the few outliers?
# (CONCLUSION: THERE IS A MUCH BETTER CORRELATION BETWEEN CONFIDENCE AND DISTANCE IN BAL)

# color the tree based on average confidence in OC errors origininating from that tip
f1 <- function(x) {
	if (is.leaf(x)) {
		m <- match(attr(x, "label"), loto$genera)
		if (is.na(m)) { # quotes in the label
			attr(x, "hits") <- 0
			attr(x, "tots") <- 0
		} else {
			attr(x, "hits") <- loto$confidences[m]
			attr(x, "edgePar") <- list(col=rgb(attr(x, "hits"), 0, 0))
			attr(x, "tots") <- 1
		}
	} else {
		x[[1]] <- f(x[[1]])
		x[[2]] <- f(x[[2]])
		
		attr(x, "hits") <- attr(x[[1]], "hits") + attr(x[[2]], "hits")
		attr(x, "tots") <- attr(x[[1]], "tots") + attr(x[[2]], "tots")
		
		attr(x, "edgePar") <- list(col=rgb(min(attr(x, "hits")/attr(x, "tots"), 1), 0, 0))
	}
	return(x)
}
loto <- loto_unbal
x <- lapply(c, f1)
attributes(x) <- attributes(c)
dev.new(height=7.75, width=38)
plot(x, nodePar=list(pch=NA, lab.cex=0.3))
loto <- loto_bal
x <- lapply(c, f1)
attributes(x) <- attributes(c)
dev.new(height=7.75, width=38)
plot(x, nodePar=list(pch=NA, lab.cex=0.3))
# (CONCLUSION: OC ERRORS ORIGINATE IN REGIONS OF TREE WITH SHORTER BRANCHES?)

# color the tree based on confidence in OC errors assigned to that tip
f2 <- function(x) {
	if (is.leaf(x)) {
		w <- which(loto$predictions==attr(x, "label"))
		if (length(w)==0) {
			attr(x, "hits") <- 0
			attr(x, "tots") <- 1
		} else {
			attr(x, "hits") <- max(loto$confidences[w])
			attr(x, "edgePar") <- list(col=rgb(attr(x, "hits"), 0, 0))
			attr(x, "tots") <- 1
		}
	} else {
		x[[1]] <- f(x[[1]])
		x[[2]] <- f(x[[2]])
		
		attr(x, "hits") <- attr(x[[1]], "hits") + attr(x[[2]], "hits")
		attr(x, "tots") <- attr(x[[1]], "tots") + attr(x[[2]], "tots")
		
		attr(x, "edgePar") <- list(col=rgb(min(attr(x, "hits")/attr(x, "tots"), 1), 0, 0))
	}
	return(x)
}
loto <- loto_unbal
x <- lapply(c, f2)
attributes(x) <- attributes(c)
dev.new(height=7.75, width=38)
plot(x, nodePar=list(pch=NA, lab.cex=0.3))
loto <- loto_bal
x <- lapply(c, f2)
attributes(x) <- attributes(c)
dev.new(height=7.75, width=38)
plot(x, nodePar=list(pch=NA, lab.cex=0.3))
# (CONCLUSION: OC ERRORS ARE ASSIGNED TO REGIONS OF TREE WITH SHORTER BRANCHES?)

