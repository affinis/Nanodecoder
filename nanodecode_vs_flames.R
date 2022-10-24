setwd('/data/LyuLin/Scripts/NanoDecoder/testcell')

flames.P.bar<-read.delim('P.flame.bcr.bar',header = F)
flames.T.tcr.bar<-read.delim('T.flame.tcr.bar',header = F)
flames.T.bcr.bar<-read.delim('T.flame.bcr.bar',header = F)
colnames(flames.P.bar)<-c("barcode_flames","UMI_flames","read")
colnames(flames.T.tcr.bar)<-c("barcode_flames","UMI_flames","read")
colnames(flames.T.bcr.bar)<-c("barcode_flames","UMI_flames","read")

nanodeco.P.bar<-read.delim('P.pro.bcr.nanodeco',header = F)
nanodeco.T.tcr.bar<-read.delim('T.pro.tcr.nanodeco',header = F)
nanodeco.T.bcr.bar<-read.delim('T.pro.bcr.nanodeco',header = F)
colnames(nanodeco.P.bar)<-c("read","flag","R1","barcode_nanodeco","barcode_pos_nanodeco","TSO","strand","mismatch")
colnames(nanodeco.T.tcr.bar)<-c("read","flag","R1","barcode_nanodeco","barcode_pos_nanodeco","TSO","strand","mismatch")
colnames(nanodeco.T.bcr.bar)<-c("read","flag","R1","barcode_nanodeco","barcode_pos_nanodeco","TSO","strand","mismatch")

P.bar<-left_join(nanodeco.P.bar,flames.P.bar,by="read")
T.bar.tcr<-left_join(nanodeco.T.tcr.bar,flames.T.tcr.bar,by="read")
T.bar.bcr<-left_join(nanodeco.T.bcr.bar,flames.T.bcr.bar,by="read")

P.bar$conflict<-ifelse(P.bar$barcode_nanodeco==P.bar$barcode_flames,"no","yes")
T.bar.tcr$conflict<-ifelse(T.bar.tcr$barcode_nanodeco==T.bar.tcr$barcode_flames,"no","yes")
T.bar.bcr$conflict<-ifelse(T.bar.bcr$barcode_nanodeco==T.bar.bcr$barcode_flames,"no","yes")

P.bar$origin<-"Plasma_BCR"
T.bar.tcr$origin<-"T_TCR"
T.bar.bcr$origin<-"T_BCR"

reads_tested<-base::Reduce(rbind,list(P.bar,T.bar.tcr,T.bar.bcr))
ggplot(reads_tested)+geom_bar(aes(x=origin,fill=flag),position="fill",color="black")+
  theme_publication()+scale_fill_d3()+scale_y_continuous(expand=c(0,0))

ggplot(reads_tested %>% dplyr::filter(.,R1>=0,R1<100))+geom_density(aes(x=R1))+
  theme_publication()+scale_y_continuous(expand=c(0,0))+xlab("R1 position")+xlim(0,100)

ggplot(reads_tested %>% dplyr::filter(.,barcode_pos_nanodeco>=0,barcode_pos_nanodeco<100))+geom_density(aes(x=barcode_pos_nanodeco))+
  theme_publication()+scale_y_continuous(expand=c(0,0))+xlab("barcode position")+xlim(0,100)

ggplot(reads_tested %>% dplyr::filter(.,TSO>=0,TSO<100))+geom_density(aes(x=TSO))+
  theme_publication()+scale_y_continuous(expand=c(0,0))+xlab("TSO position")+xlim(0,100)

ggplot(reads_tested %>% dplyr::filter(.,flag=="fine"))+geom_bar(aes(x=origin,fill=conflict),position="fill",color="black")+
  theme_publication()+scale_fill_aaas()+scale_y_continuous(expand=c(0,0))+xlab("")

ggplot(reads_tested %>% dplyr::filter(.,R1>=0,R1<100,conflict=="no"))+geom_density(aes(x=R1,color=origin))+
  theme_publication()+scale_y_continuous(expand=c(0,0))+xlab("R1 position of non-conflict reads")+xlim(0,100)

non_conflict_T_BCR<-reads_tested %>% dplyr::filter(.,R1>=0,R1<100,conflict=="no",origin=="T_BCR",flag=="fine")
non_conflict_T_BCR<-non_conflict_T_BCR$read
plotRead(Tdf,id=non_conflict_T_BCR[1])

conflict_T_BCR<-reads_tested %>% dplyr::filter(.,R1>=0,R1<100,conflict=="yes",origin=="T_BCR",flag=="fine")
conflict_T_BCR<-conflict_T_BCR$barcode_nanodeco

SCRIPTS<-"/data/LyuLin/SingleCell/SeuratLearning-2021-1-20/wrapper_script/"
source(paste(SCRIPTS,"seuratWrapper.R",sep="/"))
automatedLabeling<-function(srt,featureslist=NULL){
  selfsrt=srt
  Idents(selfsrt)="seurat_clusters"
  if(is.null(featureslist)){
    selffeatures=list("B"=c("MS4A1","BANK1","CD79A"),
                      "Stromal"=c("SPARC","CLU"),
                      "Epithelial"=c("CKB","KRT18"),
                      "Myeloid"=c("AIF1","LYZ","C1QB"),
                      "Plasma"=c("MZB1","IGHA1","JCHAIN"),
                      "TNK"=c("CD3D","GNLY","IL7R","CD2")
    )
  }else{
    selffeatures=featureslist
  }
  selfsrt=AddModuleScore(selfsrt,features=selffeatures)
  for(c in selfsrt@meta.data$seurat_clusters %>% sort %>% unique){
    df=(subset(selfsrt,idents=c))@meta.data[,paste0("Cluster",1:length(selffeatures))]
    max.cluster=colSums(df) %>% which.max()
    max.cluster.name=selffeatures[max.cluster] %>% names()
    selfsrt@meta.data[Cells(subset(selfsrt,idents=c)),"autotype"]=max.cluster.name
  }
  return(selfsrt)
}
srt<-seuratWrap1('/data/chenz/cellranger/TF2104ZJCR2X/outs/filtered_feature_bc_matrix/')
srt<-seuratWrap2(srt)
srt<-seuratWrap3(srt,res = 0.05)
DimPlot(srt)
srt<-automatedLabeling(srt)
DimPlot(srt,group.by = "autotype")
Bs<-srt@meta.data %>% dplyr::filter(.,autotype=="B") %>% rownames() %>% gsub("-1","",.)
Ps<-srt@meta.data %>% dplyr::filter(.,autotype=="Plasma") %>% rownames() %>% gsub("-1","",.)
BCR_cells<-c(Bs,Ps)
intersect(BCR_cells,conflict_T_BCR)
