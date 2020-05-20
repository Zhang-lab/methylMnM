# Rscript MnM.r MeDIP_1 MeDIP_2 MRE_1 MRE_2 output_file

ts = proc.time()

args = commandArgs(TRUE)
(MeDIP_1 = args[1])
(MeDIP_2 = args[2])
(MRE_1 = args[3])
(MRE_2 = args[4])
(Out_file = args[5])

library(methylMnM)

cpg_bin = '/expr/bzhang/Genome_related/hg19/MnM/num500_allcpg_hg19.bed'
mre_cpg_bin = '/expr/bzhang/Genome_related/hg19/MnM/num500_Five_mre_cpg_hg19.bed'
bl='/expr/bzhang/Genome_related/hg19/MnM/blacklist.bed'

Num500_MeDIP_1 = paste('num500_', MeDIP_1,sep = '')
Num500_MeDIP_2 = paste('num500_', MeDIP_2,sep = '')
Num500_MRE_1 = paste('num500_', MRE_1,sep = '')
Num500_MRE_2 = paste('num500_', MRE_2,sep = '')

countMeDIPbin(file.Medipsite = MeDIP_1,file.blacklist=bl,file.bin = cpg_bin,writefile = Num500_MeDIP_1, binlength = 500)
countMeDIPbin(file.Medipsite = MeDIP_2,file.blacklist=bl,file.bin = cpg_bin,writefile = Num500_MeDIP_2, binlength = 500)

countMREbin(file.MREsite = MRE_1,file.blacklist=bl,file.bin = cpg_bin,writefile = Num500_MRE_1, binlength = 500)
countMREbin(file.MREsite = MRE_2,file.blacklist=bl,file.bin = cpg_bin,writefile = Num500_MRE_2, binlength = 500)

# calculate p-value of each bin.

datafile = c(Num500_MeDIP_1,Num500_MeDIP_2,Num500_MRE_1,Num500_MRE_2)
pv_output = paste('pv_',Out_file,sep = '')
MnM.test(file.dataset = datafile,chrstring = NULL,file.cpgbin=cpg_bin,file.mrecpgbin=mre_cpg_bin,writefile=pv_output,mreratio=3/7,method="XXYY", psd=2,mkadded=1,a=1e-10,cut=100,top=500)

# calculate q-value of each bin
qv_output = paste('qv_',Out_file,sep = '')
MnM.qvalue(pv_output,qv_output)

qv = read.table(qv_output,header = T)
DMR_outfile = paste('DMRs_',qv_output,sep = '')
DMR_bg = paste(DMR_outfile,'.bedGraph',sep='')
dmr = MnM.selectDMR(qv,up = 1.45,down = 1/1.45,q.value = 1e-5,cutoff = 'q-value',quant = 0.9)
write.table(dmr,file=DMR_outfile,sep='\t',quote=F,row.names=F)
Dbg = dmr[,1:4]
Dbg[,4]= log10(dmr[,12]+1e-40)* (dmr[,11]/abs(dmr[,11]))
write.table(Dbg,file=DMR_bg,sep='\t',quote=F,row.names=F,col.names=F)

te = proc.time()
print('costing seconds:');te-ts
