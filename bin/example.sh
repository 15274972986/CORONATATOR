perl ../s4.annotateSG.pl -a mers.sg.anno mers.c.sg |sort -k5 -n > mers.c.sg.proc.sort
perl ../s5.clusterSG.pl mers.c.sg.proc.sort > mers.c.sg.proc.sort.cls
cat mers.c.sg.proc.sort.cls |perl ../s6.forDotplot.pl > mers.c.sg.proc.sort.cls.fdot
cat mers.c.sg | cut -f1,2|sort|uniq > mers.p.s



