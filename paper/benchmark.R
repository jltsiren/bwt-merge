# Use R --slave --args name < benchmark.R

args = commandArgs()

name = args[4]

x = 3.3
y = 3

data <- read.csv(file = paste(name, ".csv", sep = ""), head = FALSE, sep = ",", dec = ".", check.names = FALSE)
pdf(file = paste(name, ".pdf", sep = ""), width = x, height = y, paper = "special",
  family = "Helvetica", pointsize = 11)
par(mar=c(4, 4, 1, 1))

xrange = c(120, 170)
xscale = c(120, 130, 140, 150, 160, 170)
xtitle = "Memory usage (GB)"
xlabs = xscale

yrange = c(30, 36)
yscale = c(30, 31, 32, 33, 34, 35, 36)
ytitle = "Time (h)"
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

text(data[1, 3], data[1, 2] / 3600, data[1, 1], cex = 0.8, pos = 2)
text(data[2, 3], data[2, 2] / 3600, data[2, 1], cex = 0.8, pos = 1)
text(data[3, 3], data[3, 2] / 3600, data[3, 1], cex = 0.8, pos = 1)
text(data[4, 3], data[4, 2] / 3600, data[4, 1], cex = 0.8, pos = 2)
text(data[5, 3], data[5, 2] / 3600, data[5, 1], cex = 0.8, adj = c(0.65, 2))
text(data[6, 3], data[6, 2] / 3600, data[6, 1], cex = 0.8, pos = 2)
text(data[7, 3], data[7, 2] / 3600, data[7, 1], cex = 0.8, pos = 2)
text(data[8, 3], data[8, 2] / 3600, data[8, 1], cex = 0.8, pos = 2)

segments(124.2, 29, 124.2, 37, lty = "dashed")

dev.off()
q()
