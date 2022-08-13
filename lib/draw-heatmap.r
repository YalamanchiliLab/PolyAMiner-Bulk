library(gplots)

# Read normalized counts from the standard input.
data = read.table("stdin", header = T, sep = "\t", as.is = TRUE)

# The gene names in the first column.
gene = data[, 1]

# Load the data from the second column on.
vals = as.matrix(data[, 2 : ncol(data)])

# Adds a little noise to each element to avoid the
# clustering function fail on zero variance datalines.
vals = jitter(vals, factor = 1, amount = 0.00001)

# Each row is normalized to a z-score
zscore = NULL
for (i in 1 : nrow(vals)) {
    row = vals[i,]
    zrow = (row - mean(row)) / sd(row)
    zscore = rbind(zscore, zrow)
}

# Add back gene names as row names.
row.names(zscore) = gene

# Turn it into a matrix for heatmap2.
zscore = as.matrix(zscore)

# Open the drawing device.
pdf('|cat')

# Set the color scheme.
colors = colorRampPalette(c("green", "black", "red"), space = "rgb")(256)

# Draw the heatmap.
heatmap.2(zscore, col = colors, density.info = "none", trace = "none", margins = c(16, 7), lhei = c(1, 5), cexCol = 1, dendrogram = "column", key = FALSE, labRow = FALSE)

# Turn off the device.
dev.off()

