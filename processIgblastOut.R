library("phangorn")
library("ggtree")
library("treeio")
library("ggimage")
library("TDbook")
library("seqinr")
library("Biostrings")
library("ggseqlogo")

sample_prefix<-"test"
immune_molecule_prefix<-"Ig"
visium_out_path<-"/data/LyuLin/Spatial/fastqs/nano/M1054/ST2201CLLX4X"

Igblast_out_Ig<-read.delim('/data/LyuLin/Spatial/fastqs/nano/M1054/TRUST4/igblast_out/Igs.tsv')
Igblast_out_TR<-read.delim('/data/LyuLin/Spatial/fastqs/nano/M1054/TRUST4/igblast_out/TRs.tsv')
ggplot(Igblast_out_Ig)+geom_violin(aes(x=locus,y=v_identity))+scale_x_discrete(labels = c("unknown", "IGH", "IGK", "IGL"))+
  theme_publication(base.size = 16)+scale_y_continuous(expand=c(0,0))
ggplot(Igblast_out_TR)+geom_violin(aes(x=locus,y=v_identity))+scale_x_discrete(labels = c("unknown", "TRA", "TRB", "TRD","TRG"))+
  theme_publication(base.size = 16)+scale_y_continuous(expand=c(0,0))
idents<-c(60,65,70,75,80,85,90,95)
shared_reads<-NULL
uniq_ig_reads<-NULL
uniq_tcr_reads<-NULL
for (ident in idents) {
  ig_reads=dplyr::filter(Igblast_out_Ig,!is.na(v_identity),!is.na(j_identity),v_identity>ident,j_identity>ident)$sequence_id
  tcr_reads=dplyr::filter(Igblast_out_TR,!is.na(v_identity),!is.na(j_identity),v_identity>ident,j_identity>ident)$sequence_id
  shared_read=intersect(ig_reads,tcr_reads) %>% length()
  shared_reads=append(shared_reads,shared_read)
  
  uniq_ig_read=setdiff(ig_reads,tcr_reads) %>% length()
  uniq_tcr_read=setdiff(tcr_reads,ig_reads) %>% length()
  uniq_ig_reads=append(uniq_ig_reads,uniq_ig_read)
  uniq_tcr_reads=append(uniq_tcr_reads,uniq_tcr_read)
}
identity_vs_shared<-data.frame(identity_set=idents,num_of_shared_reads=shared_reads,uniquely_mapped_to_Ig=uniq_ig_reads,uniquely_mapped_to_TCR=uniq_tcr_reads)
identity_vs_shared<-gather(identity_vs_shared,key="readstat",value="number_of_reads",-identity_set)
ggplot(identity_vs_shared)+geom_path(aes(x=identity_set,y=number_of_reads,group=readstat,color=readstat))+
  theme_publication(base.size = 16)+scale_color_aaas()

Igblast_out_Ig<-dplyr::filter(Igblast_out_Ig,!is.na(v_identity),!is.na(j_identity),v_identity>85,j_identity>85)
Igblast_out_Ig[is.na(Igblast_out_Ig)]<-"unknown/NA"
ggplot(Igblast_out_Ig)+geom_bar(aes(x=locus),fill=pal_d3()(10)[1])+
  theme_publication(base.size = 16)+scale_y_continuous(expand=c(0,0))

Igblast_out_TR<-dplyr::filter(Igblast_out_TR,!is.na(v_identity),!is.na(j_identity),v_identity>85,j_identity>85)
Igblast_out_TR[is.na(Igblast_out_TR)]<-"unknown/NA"
ggplot(Igblast_out_TR)+geom_bar(aes(x=locus),fill=pal_d3()(10)[2])+
  theme_publication(base.size = 16)+scale_y_continuous(expand=c(0,0))

Igblast_out_Ig$note<-"Ig_seq"
Igblast_out_TR$note<-"TCR_seq"

ggplot(Igblast_out_Ig)+geom_bar(aes(x=note,fill=locus),position="fill")+xlab("")+
  theme_publication(base.size = 16)+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=pal_d3()(10)[3:5])
ggplot(Igblast_out_TR)+geom_bar(aes(x=note,fill=locus),position="fill")+xlab("")+
  theme_publication(base.size = 16)+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=pal_d3()(10)[6:9])

Igblast_summary<-read.delim('/data/LyuLin/Spatial/fastqs/nano/M1054/TRUST4/igblast_out/Ig_clone_summary.tsv',header = F)
Igblast_summary<-Igblast_summary %>% dplyr::filter(.,V7=="Yes")
productive_counts<-ggplot(Igblast_summary)+geom_bar(aes(x=fct_inorder(V1),y=V3),stat="identity")+theme_publication(base.size=16)+
  scale_y_continuous(expand=c(0,0))+theme(axis.text.y=element_blank(),axis.ticks.y = element_blank())+
  ylab("count")+xlab("")+coord_flip()

# read the MSA data
msa<-read.fasta("/data/LyuLin/Spatial/fastqs/nano/M1054/TRUST4/igblast_out/cdr3_MSA/cdr3_proliferable.aln.fasta")
msa_names<-names(msa)
msa_names<-msa_names %>% gsub("*","_",.,fixed = T)
msa.df<-lapply(msa,as.vector) %>% as.data.frame() 
colnames(msa.df)<-msa_names
msa.df$position<-rownames(msa.df)
msa.df<-gather(msa.df,key="cloneType",value="base",-position)
msa.df$base<-toupper(msa.df$base)
productive_base_usage<-ggplot(msa.df)+geom_tile(aes(x=fct_inorder(cloneType),y=fct_inorder(position),fill=base))+ylab("position in CDR3")+
  scale_fill_jama()+theme_publication(base.size=16)+xlab("clone_type")+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
    axis.text.x=element_blank(),axis.ticks.x=element_blank())+coord_flip()
get_legend(productive_base_usage) %>% plot()
productive_base_usage<-productive_base_usage+NoLegend()

msa<-read.fasta("/data/LyuLin/Spatial/fastqs/nano/M1054/TRUST4/igblast_out/cdr3_MSA/cdr3_proliferable.aln.faa")
msa_names<-names(msa)
msa_names<-msa_names %>% gsub("*","_",.,fixed = T)
msa.aa.df<-lapply(msa,as.vector) %>% as.data.frame()
colnames(msa.aa.df)<-msa_names
msa.aa.df$position<-rownames(msa.aa.df)
msa.aa.df<-gather(msa.aa.df,key="cloneType",value="AA",-position)
msa.aa.df$AA<-toupper(msa.aa.df$AA)
productive_aa_usage<-ggplot(msa.aa.df)+geom_tile(aes(x=fct_inorder(cloneType),y=fct_inorder(position),fill=AA))+
  ylab("position in CDR3")+theme_publication(base.size=16)+xlab("")+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  coord_flip()+scale_fill_manual(values=c(pal_d3("category20b")(20),pal_d3()(10)))
get_legend(productive_aa_usage) %>% plot()
productive_aa_usage<-productive_aa_usage+NoLegend()

ggarrange(productive_base_usage,productive_aa_usage,productive_counts,align = "h",ncol=3,common.legend = F,widths = c(2,1,1))
ggsave("/data/LyuLin/Spatial/fastqs/nano/M1054/TRUST4/igblast_out/cdr3_MSA/base_usage.pdf",width=12,height=12)


##generate Ig srt
Igblast_out_Ig$read_id<-strsplit(Igblast_out_Ig$sequence_id,"_") %>% lapply(.,`[`,1) %>% unlist()
Igblast_out_Ig$barcode<-strsplit(Igblast_out_Ig$sequence_id,"_") %>% lapply(.,`[`,2) %>% unlist()
Igblast_out_Ig<-Igblast_out_Ig[,c("read_id","barcode")]
rownames(Igblast_out_Ig)<-Igblast_out_Ig$read_id

Igblast_clones<-read.delim('/data/LyuLin/Spatial/fastqs/nano/M1054/TRUST4/igblast_out/Ig_clone2reads',header = F)
Igblast_clones<-Igblast_clones[,c(1,7)]
rownames(Igblast_clones)<-Igblast_clones$V1
colnames(Igblast_clones)<-c("clone_id","reads")
clone_summary<-data.frame("barcode"=as.character(NULL),"count"=as.numeric(NULL),"clone_id"=as.character(NULL))
for(rown in rownames(Igblast_clones)){
  reads=Igblast_clones[rown,"reads"]
  reads_id=strsplit(reads,",") %>% unlist() %>% strsplit(.,"_") %>% lapply(.,`[`,1) %>% unlist()
  counts=Igblast_out_Ig[intersect(reads_id,rownames(Igblast_out_Ig)),] %>% group_by(barcode) %>% summarise(.,count=n())
  counts$clone_id=paste0(immune_molecule_prefix,"_clone",rown)
  clone_summary=rbind(clone_summary,counts)
}
clone_summary$barcode<-paste0(clone_summary$barcode,"-1")
clone_mtx<-spread(clone_summary,key=clone_id,value=count)
clone_mtx<-as.data.frame(clone_mtx)
rownames(clone_mtx)<-clone_mtx$barcode
clone_mtx$barcode<-NULL
clone_mtx<-t(clone_mtx)
clone_mtx[is.na(clone_mtx)]<-0
srt<-Load10X_Spatial(visium_out_path,filter.matrix=F,filename="raw_feature_bc_matrix.h5")
allbarcodes=Cells(srt)
invalidbarcodes=allbarcodes[!(allbarcodes %in% colnames(clone_mtx))]
zeros=data.frame(row.names = row.names(clone_mtx),matrix(0,nrow = nrow(clone_mtx),ncol = length(invalidbarcodes)))
colnames(zeros)=invalidbarcodes
clone_mtx=cbind(clone_mtx,zeros)
clone_mtx=clone_mtx[,allbarcodes]
assay_SIP=CreateAssayObject(clone_mtx)
srt@assays$SIP=assay_SIP
DefaultAssay(srt)="SIP"
srt@assays$SIP@key="spatial_"
srt=AddMetaData(srt,metadata=colSums(srt),col.name="nCount_Clones")
nfeature=colSums(srt@assays$SIP@counts>0)
srt=AddMetaData(srt,metadata=nfeature,col.name="nFeature_Clones")

SpatialFeaturePlot(srt,"nFeature_Clones")
saveRDS(srt,paste0("test_Ig_srt.rds"))


##generate TCR srt
Igblast_out_TR$read_id<-strsplit(Igblast_out_TR$sequence_id,"_") %>% lapply(.,`[`,1) %>% unlist()
Igblast_out_TR$barcode<-strsplit(Igblast_out_TR$sequence_id,"_") %>% lapply(.,`[`,2) %>% unlist()
Igblast_out_TR<-Igblast_out_TR[,c("read_id","barcode")]
rownames(Igblast_out_TR)<-Igblast_out_TR$read_id

Igblast_clones<-read.delim('/data/LyuLin/Spatial/fastqs/nano/M1054/TRUST4/igblast_out/TCR_clone2reads',header = F)
Igblast_clones<-Igblast_clones[,c(1,7)]
rownames(Igblast_clones)<-Igblast_clones$V1
colnames(Igblast_clones)<-c("clone_id","reads")
clone_summary<-data.frame("barcode"=as.character(NULL),"count"=as.numeric(NULL),"clone_id"=as.character(NULL))
for(rown in rownames(Igblast_clones)){
  reads=Igblast_clones[rown,"reads"]
  reads_id=strsplit(reads,",") %>% unlist() %>% strsplit(.,"_") %>% lapply(.,`[`,1) %>% unlist()
  counts=Igblast_out_TR[intersect(reads_id,rownames(Igblast_out_TR)),] %>% group_by(barcode) %>% summarise(.,count=n())
  counts$clone_id=paste0(immune_molecule_prefix,"_clone",rown)
  clone_summary=rbind(clone_summary,counts)
}
clone_summary$barcode<-paste0(clone_summary$barcode,"-1")
clone_mtx<-spread(clone_summary,key=clone_id,value=count)
clone_mtx<-as.data.frame(clone_mtx)
rownames(clone_mtx)<-clone_mtx$barcode
clone_mtx$barcode<-NULL
clone_mtx<-t(clone_mtx)
clone_mtx[is.na(clone_mtx)]<-0
srt<-Load10X_Spatial(visium_out_path,filter.matrix=F,filename="raw_feature_bc_matrix.h5")
allbarcodes=Cells(srt)
invalidbarcodes=allbarcodes[!(allbarcodes %in% colnames(clone_mtx))]
zeros=data.frame(row.names = row.names(clone_mtx),matrix(0,nrow = nrow(clone_mtx),ncol = length(invalidbarcodes)))
colnames(zeros)=invalidbarcodes
clone_mtx=cbind(clone_mtx,zeros)
clone_mtx=clone_mtx[,allbarcodes]
assay_SIP=CreateAssayObject(clone_mtx)
srt@assays$SIP=assay_SIP
DefaultAssay(srt)="SIP"
srt@assays$SIP@key="spatial_"
srt=AddMetaData(srt,metadata=colSums(srt),col.name="nCount_Clones")
nfeature=colSums(srt@assays$SIP@counts>0)
srt=AddMetaData(srt,metadata=nfeature,col.name="nFeature_Clones")

SpatialFeaturePlot(srt,"nCount_Clones")
saveRDS(srt,paste0("test_TCR_srt.rds"))

if(!("srt" %in% ls())){
  srt<-readRDS("/data/LyuLin/Spatial/fastqs/nano/M1054/TRUST4/igblast_out/test_Ig_srt.rds")
}
if(!("immune_molecule_prefix" %in% ls())){
  immune_molecule_prefix<-"Ig"
}
Igblast_summary<-read.delim('/data/LyuLin/Spatial/fastqs/nano/M1054/TRUST4/igblast_out/Ig_clone_summary.tsv',header = F)
productive_type<-Igblast_summary[Igblast_summary$V7=="Yes","V1"]
productive_type<-productive_type %>% paste0("Ig-clone",.)
noproductive_type<-Igblast_summary[Igblast_summary$V7=="No","V1"]
noproductive_type<-noproductive_type %>% paste0("Ig-clone",.)
nCount_productive<-FetchData(srt,vars=productive_type,cells=Cells(srt)) %>% rowSums()
nCount_nonproductive<-FetchData(srt,vars=noproductive_type,cells=Cells(srt)) %>% rowSums()
srt<-AddMetaData(srt,nCount_productive,col.name="nCount_Productive")
srt<-AddMetaData(srt,nCount_nonproductive,col.name="nCount_nonProductive")
SpatialFeaturePlot(srt,"nCount_Productive")
SpatialFeaturePlot(srt,"nCount_nonProductive")



#######################
#using igblast raw out
Igblast_out_TR<-read.delim('/data/LyuLin/Spatial/fastqs/nano/M1054/TRUST4/igblast_out/TRs.tsv')
Igblast_out_TR$read_id<-Igblast_out_TR$sequence_id %>% strsplit(.,"_") %>% lapply(.,`[`,1) %>% unlist()
Igblast_out_TR$barcode<-Igblast_out_TR$sequence_id %>% strsplit(.,"_") %>% lapply(.,`[`,2) %>% unlist()
rownames(Igblast_out_TR)<-Igblast_out_TR$read_id
Igblast_out_TR<-Igblast_out_TR[,c("read_id","barcode")]
TCR_clone2reads<-read.delim('/data/LyuLin/Spatial/fastqs/nano/M1054/TRUST4/igblast_out/TCR_clone2reads',header = F)
clone_id_tcr<-TCR_clone2reads$V1
TCR_id_bar_clone<-data.frame("read_id"=as.character(NULL),"barcode"=as.numeric(NULL),"clone_id"=as.character(NULL))
rownames(TCR_clone2reads)<-TCR_clone2reads$V1
for(clone in clone_id_tcr){
  reads=TCR_clone2reads[clone,"V7"]
  reads=strsplit(reads,",") %>% unlist() %>% strsplit(.,"_") %>% lapply(.,`[`,1) %>% unlist()
  df=Igblast_out_TR[reads,]
  df$clone_id=clone
  TCR_id_bar_clone=rbind(TCR_id_bar_clone,df)
}
TCR_id_bar_clone=TCR_id_bar_clone %>% group_by(barcode,clone_id) %>% summarise(.,count=n())
TCR_id_bar_clone$barcode<-paste0(TCR_id_bar_clone$barcode,"-1")
clone_mtx<-spread(TCR_id_bar_clone,key=clone_id,value=count)
clone_mtx<-as.data.frame(clone_mtx)
rownames(clone_mtx)<-clone_mtx$barcode
clone_mtx$barcode<-NULL
clone_mtx<-t(clone_mtx)
clone_mtx[is.na(clone_mtx)]<-0
srt<-Load10X_Spatial(visium_out_path,filter.matrix=F,filename="raw_feature_bc_matrix.h5")
allbarcodes=Cells(srt)
invalidbarcodes=allbarcodes[!(allbarcodes %in% colnames(clone_mtx))]
zeros=data.frame(row.names = row.names(clone_mtx),matrix(0,nrow = nrow(clone_mtx),ncol = length(invalidbarcodes)))
colnames(zeros)=invalidbarcodes
clone_mtx=cbind(clone_mtx,zeros)
clone_mtx=clone_mtx[,allbarcodes]
assay_SIP=CreateAssayObject(clone_mtx)
srt@assays$SIP=assay_SIP
DefaultAssay(srt)="SIP"
srt@assays$SIP@key="spatial_"
srt=AddMetaData(srt,metadata=colSums(srt),col.name="nCount_Clones")
nfeature=colSums(srt@assays$SIP@counts>0)
srt=AddMetaData(srt,metadata=nfeature,col.name="nFeature_Clones")

TRs.valid<-read.delim('/data/LyuLin/Spatial/fastqs/nano/M1054/TRUST4/igblast_out/TRs.valid.tsv')
TRs.valid<-TRs.valid[,c(1,3)]
colnames(TRs.valid)<-c("read_id","locus")
ggplot(TRs.valid)+geom_bar(aes(x=locus),fill=pal_d3()(10)[2])+
  theme_publication(base.size = 16)+scale_y_continuous(expand=c(0,0))
TRs.valid$note<-"TCR_seq"
ggplot(TRs.valid)+geom_bar(aes(x=note,fill=locus),position="fill")+xlab("")+
  theme_publication(base.size = 16)+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=pal_d3()(10)[6:9])

Igblast_out_Ig<-read.delim('/data/LyuLin/Spatial/fastqs/nano/M1054/TRUST4/igblast_out/Igs.tsv')
Igblast_out_Ig$read_id<-Igblast_out_Ig$sequence_id %>% strsplit(.,"_") %>% lapply(.,`[`,1) %>% unlist()
Igblast_out_Ig$barcode<-Igblast_out_Ig$sequence_id %>% strsplit(.,"_") %>% lapply(.,`[`,2) %>% unlist()
rownames(Igblast_out_Ig)<-Igblast_out_Ig$read_id
Igblast_out_Ig<-Igblast_out_Ig[,c("read_id","barcode")]
Ig_clone2reads<-read.delim('/data/LyuLin/Spatial/fastqs/nano/M1054/TRUST4/igblast_out/Ig_clone2reads',header = F)
clone_id_ig<-Ig_clone2reads$V1
IG_id_bar_clone<-data.frame("read_id"=as.character(NULL),"barcode"=as.numeric(NULL),"clone_id"=as.character(NULL))
rownames(Ig_clone2reads)<-Ig_clone2reads$V1
for(clone in clone_id_ig){
  reads=Ig_clone2reads[clone,"V7"]
  reads=strsplit(reads,",") %>% unlist() %>% strsplit(.,"_") %>% lapply(.,`[`,1) %>% unlist()
  df=Igblast_out_Ig[reads,]
  df$clone_id=clone
  IG_id_bar_clone=rbind(IG_id_bar_clone,df)
}
IG_id_bar_clone=IG_id_bar_clone %>% group_by(barcode,clone_id) %>% summarise(.,count=n())
IG_id_bar_clone$barcode<-paste0(IG_id_bar_clone$barcode,"-1")
clone_mtx<-spread(IG_id_bar_clone,key=clone_id,value=count)
clone_mtx<-as.data.frame(clone_mtx)
rownames(clone_mtx)<-clone_mtx$barcode
clone_mtx$barcode<-NULL
clone_mtx<-t(clone_mtx)
clone_mtx[is.na(clone_mtx)]<-0
srt2<-Load10X_Spatial(visium_out_path,filter.matrix=F,filename="raw_feature_bc_matrix.h5")
allbarcodes=Cells(srt2)
invalidbarcodes=allbarcodes[!(allbarcodes %in% colnames(clone_mtx))]
zeros=data.frame(row.names = row.names(clone_mtx),matrix(0,nrow = nrow(clone_mtx),ncol = length(invalidbarcodes)))
colnames(zeros)=invalidbarcodes
clone_mtx=cbind(clone_mtx,zeros)
clone_mtx=clone_mtx[,allbarcodes]
assay_SIP=CreateAssayObject(clone_mtx)
srt2@assays$SIP=assay_SIP
DefaultAssay(srt2)="SIP"
srt2@assays$SIP@key="spatial_"
srt2=AddMetaData(srt2,metadata=colSums(srt2),col.name="nCount_Clones")
nfeature=colSums(srt2@assays$SIP@counts>0)
srt2=AddMetaData(srt2,metadata=nfeature,col.name="nFeature_Clones")

Igs.valid<-read.delim('/data/LyuLin/Spatial/fastqs/nano/M1054/TRUST4/igblast_out/Igs.valid.tsv')
Igs.valid<-Igs.valid[,c(1,3)]
colnames(Igs.valid)<-c("read_id","locus")
ggplot(Igs.valid)+geom_bar(aes(x=locus),fill=pal_d3()(10)[1])+
  theme_publication(base.size = 16)+scale_y_continuous(expand=c(0,0))
Igs.valid$note<-"IG_seq"
ggplot(Igs.valid)+geom_bar(aes(x=note,fill=locus),position="fill")+xlab("")+
  theme_publication(base.size = 16)+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=pal_d3()(10)[3:5])

saveRDS(srt,"/data/LyuLin/Spatial/fastqs/nano/M1054/TRUST4/igblast_out/igblastRaw_TCR_srt.rds")
saveRDS(srt2,"/data/LyuLin/Spatial/fastqs/nano/M1054/TRUST4/igblast_out/igblastRaw_Ig_srt.rds")

