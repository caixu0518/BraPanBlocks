library(ggplot2)

mytheme2 <- theme_bw() + theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),

                               axis.title.x = element_blank(),
                               axis.title.y = element_text(size=8,angle=90),
                               #axis.title.y = element_text(size=8,angle=90),
                               axis.text.x=element_text(size=9,angle=0),
                               axis.text.y = element_text(size=7),  ##- biao qian
                               #axis.ticks.y = element_blank(), ##- ke du xian 
                               legend.position="none",
                               legend.key.size = unit(0.4,'cm'),
                               legend.text = element_text(size=6,angle=0),
                               strip.text = element_text(size=8)

                              )

argv<-commandArgs(TRUE)

data1 <- read.table(argv[1], header = T, sep="\t")

data1$Chr <- factor(data1$Chr, levels=c("A01", "A02", "A03", "A04", "A05", "A06", "A07", "A08", "A09", "A10"))

chrColor <- c("#FFFFFF", "#9933FF", "#FF33FF", "#FF3399", "#FF3333", "#FF9933", "#99FF33", "#33FF99", "#3399FF", "#B8B800")

p2 <- ggplot(data1, aes(x=Pos/1000000, y=zFst, color=Chr, fill=Chr)) +
      geom_point(size=0.3, shape=1, alpha=0.5) +
      ylim(0, 7) +
      facet_grid(. ~ Chr,as.table=TRUE, scales="free_x", space="free_x")  +
       
      #geom_vline(data=subset(data1, Chr=="A01"), aes(xintercept=16000001/1000000), colour="#990000", linetype="dashed")  +
      #geom_vline(data=subset(data1, Chr=="A01"), aes(xintercept=23900000/1000000), colour="#990000", linetype="dashed")  +

      #geom_vline(data=subset(data1, Chr=="A02"), aes(xintercept=18600001/1000000), colour="#990000", linetype="dashed")  +
      #geom_vline(data=subset(data1, Chr=="A02"), aes(xintercept=21500000/1000000), colour="#990000", linetype="dashed")  +

      #geom_vline(data=subset(data1, Chr=="A03"), aes(xintercept=34400001/1000000), colour="#990000", linetype="dashed")  +
      #geom_vline(data=subset(data1, Chr=="A03"), aes(xintercept=39500000/1000000), colour="#990000", linetype="dashed")  +

      #geom_vline(data=subset(data1, Chr=="A04"), aes(xintercept=6500001/1000000), colour="#990000", linetype="dashed")  +
      #geom_vline(data=subset(data1, Chr=="A04"), aes(xintercept=6900000/1000000), colour="#990000", linetype="dashed")  +

      #geom_vline(data=subset(data1, Chr=="A05"), aes(xintercept=13300001/1000000), colour="#990000", linetype="dashed")  +
      #geom_vline(data=subset(data1, Chr=="A05"), aes(xintercept=29400000/1000000), colour="#990000", linetype="dashed")  +

      #geom_vline(data=subset(data1, Chr=="A06"), aes(xintercept=15000001/1000000), colour="#990000", linetype="dashed")  +
      #geom_vline(data=subset(data1, Chr=="A06"), aes(xintercept=45800000/1000000), colour="#990000", linetype="dashed")  +

      #geom_vline(data=subset(data1, Chr=="A07"), aes(xintercept=5300001/1000000), colour="#990000", linetype="dashed")  +
      #geom_vline(data=subset(data1, Chr=="A07"), aes(xintercept=6800000/1000000), colour="#990000", linetype="dashed")  +

      #geom_vline(data=subset(data1, Chr=="A08"), aes(xintercept=5900001/1000000), colour="#990000", linetype="dashed")  +
      #geom_vline(data=subset(data1, Chr=="A08"), aes(xintercept=14800000/1000000), colour="#990000", linetype="dashed")  +

      #geom_vline(data=subset(data1, Chr=="A09"), aes(xintercept=21700001/1000000), colour="#990000", linetype="dashed")  +
      #geom_vline(data=subset(data1, Chr=="A09"), aes(xintercept=46600000/1000000), colour="#990000", linetype="dashed")  +

      #geom_vline(data=subset(data1, Chr=="A10"), aes(xintercept=4000001/1000000), colour="#990000", linetype="dashed")  +
      #geom_vline(data=subset(data1, Chr=="A10"), aes(xintercept=7100000/1000000), colour="#990000", linetype="dashed")  +

      geom_hline(aes(yintercept=2.33), colour="#990000", linetype="dashed")  + 
     
      ylab("zFst") + mytheme2

ggsave(filename=argv[2],plot=p2,height=4,width=18)



