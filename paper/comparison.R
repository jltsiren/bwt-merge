# Use R --slave --args name < comparison.R

args = commandArgs()

name = args[4]

x = 3
y = 3

data <- read.csv(file = paste(name, ".csv", sep = ""), head = FALSE, sep = ",", dec = ".", check.names = FALSE)
pdf(file = paste(name, ".pdf", sep = ""), width = x, height = y, paper = "special",
  family = "Helvetica", pointsize = 11)
par(mar=c(4, 4, 1, 1))

xrange = c(0, 256)
xscale = c(0, 64, 128, 192, 256)
xtitle = "Memory usage (GB)"
xlabs = xscale

yrange = c(0, 72)
yscale = c(0, 12, 24, 36, 48, 60, 72)
ytitle = ""
ylabs = yscale

plot(c(1),
  c(1),
  type = "n",
  axes = F,
  main = "",
  xlab = xtitle,
  ylab = ytitle,
  xlim = xrange,
  ylim = yrange)

axis(1, at = xscale, lab = xlabs, cex.axis = 0.8)
axis(2, at = yscale, lab = ylabs, cex.axis = 0.8)
box()

nr = nrow(data)
nc = ncol(data)

points(data[1:nr, 3], data[1:nr, 2] / 3600, type = "p", pch = 20)

text(data[1, 3], data[1, 2] / 3600, data[1, 1], cex = 0.8, pos = 2) # RopeBWT
text(data[2, 3], data[2, 2] / 3600, data[2, 1], cex = 0.8, pos = 2) # RopeBWT2
text(data[3, 3], data[3, 2] / 3600, data[3, 1], cex = 0.8, pos = 3) # BWT-merge
text(data[4, 3], data[4, 2] / 3600, data[4, 1], cex = 0.8, pos = 1) # RopeBWT (RLO)
text(data[5, 3], data[5, 2] / 3600, data[5, 1], cex = 0.8, pos = 1) # RopeBWT2 (RLO)
text(data[6, 3], data[6, 2] / 3600, data[6, 1], cex = 0.8, pos = 1) # BWT-merge (RLO)

dev.off()
q()
