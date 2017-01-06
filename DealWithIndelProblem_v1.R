k <- 2 # k-mer size

train <- "ACTTGAACCCCCCGGAC"
test <- "TTGAACCGG"

# split into k-mers
trainKmers <- substring(train, 1:(nchar(train) - k + 1), k:nchar(train))
testKmers <- substring(test, 1:(nchar(test) - k + 1), k:nchar(test))

# replace duplicates with NAs to prevent matching
trainKmers[trainKmers %in% trainKmers[duplicated(trainKmers)]] <- NA_character_
testKmers[testKmers %in% testKmers[duplicated(testKmers)]] <- NA_character_

# find the first and last match
m <- rep(NA_integer_, length(testKmers))
m[!is.na(testKmers)] <- match(testKmers[!is.na(testKmers)], trainKmers)

# Option A:  eliminate pairs of k-mer matches that are more separated in the training sequence

eliminate <- logical(length(m))
w <- which(!is.na(m))
for (i in seq_len(length(w) - 1)) {
	if ((m[w[i + 1]] - m[w[i]]) > (w[i + 1] - w[i])) {
		eliminate[w[i + 1]] <- TRUE
		eliminate[w[i]] <- TRUE
	}
}
m[eliminate] <- NA_integer_

print(eliminate)
matches <- which(!is.na(m) & !eliminate)
print(length(matches)) # the number of hits after indel correction

# Option B:  find the boundaries of the test sequence in the training sequence


start <- NA_integer_
for (i in seq_along(m)) {
	if (!is.na(m[i])) {
		start <- m[i] - i + 1L
		break
	}
}
m <- rev(m)
end <- NA_integer_
for (i in seq_along(m)) {
	if (!is.na(m[i])) {
		end <- m[i] + i + k - 2L
		break
	}
}
if (is.na(start)) {
	cat("No matches! Don't bother testing this pair.")
} else {
	if (start < 1L)
		start <- 1L
	if (end > nchar(train))
		end <- nchar(train)
	cat("Start in train =", start, "End in train =", end)
}
substring(train, start, end)


