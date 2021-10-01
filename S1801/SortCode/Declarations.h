TList *tigList, *s3List, *tigtigList, *tigs3List, *emmaList, *icList, *tigemmaList, *tgList, *bgList, *PIDgatedList, *PIDgatedListS3; 

//Raw TIGRESS
TH1F *tigE, *addE;
TH2F *tigE_ANum,*addE_ANum, *num_addr; 

//Raw S3
TH1F *s3_E, *s3_rings_sectors_singlesT, *s3_ringNum, *s3_secNum;
TH2F *rings, *sectors, *s3_E_theta, *hitmap, *s3_rings_sectors, *s3_rings_sectors_singles, *s3_rings_sectors_gated, *s3_rings_sectors_singlesTvT, *s3_rings_sectors_singlesTvSec, *s3_rings_sectors_singlesTvRing, *s3_rings_sectors_singlesTvE;

//TIGRESS-TIGRESS
TH1F *addT_addT; 
TH2F *addE_addE, *addE_addE_tg;

//TIGRESS-S3 
TH1F *addT_s3T, *excE, *addDopp, *s3_E_gated;
TH2F *excE_theta, *addDopp_ANum, *addDopp_exc, *s3_E_theta_gated, *addE_s3_E, *addDopp_s3_E; 

//Raw EMMA
TH1F *xPos, *yPos, *siE, *icSum, *icE[5], *ssbE[2];
TH2F *pgac, *icN, *iC0V1, *iC0V2, *iC0V3, *iC1V2, *iC1V3, *iC2V3, *tDE, *siET, *ssbET[2],*icSumVSi, *ic0VSi, *ic1VSi, *ic2VSi, *ic3VSi; 

//ToF Spectra
TH1F *excE_tg, *addDopp_tg, *excE_bg, *addDopp_bg, *addemmatof, *s3emmatof;
TH2F *addE_tof, *excE_theta_tg, *addDopp_ANum_tg, *addE_addE_tofg, *addDopp_exc_tg, *addDopp_addDopp_tg, *excE_theta_bg, *addDopp_ANum_bg, *addDopp_exc_bg, *hitmap_time_gated, *icSumVSi_gated_350; 

//Mass Spectra
TH1F *excE_PIDG, *addDopp_PIDG, *excE_PIDG_bg, *addDopp_PIDG_bg, *aE_icSumVSi_gated;
TH2F *excE_theta_PIDG, *addDopp_ANum_PIDG, *addDopp_exc_PIDG, *excE_theta_PIDG_bg, *addDopp_ANum_PIDG_bg, *addDopp_exc_PIDG_bg, *rings_PIDG, *s3_E_theta_gated_PIDG; 

//TCuts
TCutG *massGate, *siicGate;
