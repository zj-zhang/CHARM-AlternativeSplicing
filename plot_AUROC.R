library("PRROC")

args = commandArgs(T)
CLFDIR = args[1]
pdf(file.path(CLFDIR, "02-test_auroc.pdf"), width=5, height=5)
data = read.table(file.path(CLFDIR, "duke_pred.data.txt"), sep="\t", header=T)
duke_pred = data$duke_pred
duke_y = data$duke_y
roc = roc.curve(as.vector(duke_pred[duke_y==1]), as.vector(duke_pred[duke_y==0]), curve=T)
plot(roc, main='Duke prediction AUROC')
dev.off()
