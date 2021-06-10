#define NTIP 128 //number of TIP channels
#define NTIPRING 10 //number of TIP rings
#define NTIGRING 6 //number of TIGRESS rings

TList *tigList, *tipList, *tipPIDList, *tipPIDGateList, *tipPIDGatedList, *tigtigList, *tigbgoList, *tiptigList; 

//Raw TIGRESS
TH1F *tigE, *tigE_unsupp, *addE, *addE_ring[NTIGRING];
TH2F *tigE_ANum,*addE_ANum, *addE_theta, *num_addr; 
TH1F *addE_ANum_Downstream,*addE_ANum_Upstream, *addE_ANum_Corona;
TH2I *tigTigDetectors, *tigTigCloverDetectors;

//TIP
TH1F *tip_E, *tip_CFDFitDiff, *tip_CFDFitDiffHighE, *tip_wfrmsize, *tiptipT;
TH2F *tip_E_pos, *tip_E_waveformAmp, *tip_lastEvtTDiff_pos;
TH1I *tip_mult, *tip_pos, *tip_ring;
TH2I *tipTipDetectors;
TH1F *tipPileupVsRing;

//TIP PID
TH2F *tip_E_PID_Sum, *tip_E_PID[NTIP];

//TIP PID Gated
TH2I *tip_particle_mult;

//TIGRESS-TIGRESS
TH1F *addT_addT;
TH2F *addE_addE, *addE_addE_tg;

//TIGRESS-BGO
TH1I *bgo_mult, *bgo_det;
TH1F *tigT_bgoT, *tigT_bgoT_supp;

//TIGRESS-TIP
TH2I *tiptig_mult;
TH1F *tipT_tigT_diff, *tipTCFD_tigT_diff, *tigE_TIPtg, *tipE_tigtg, *tipT_tigT_difftg, *tipTCFD_tigT_difftg;
TH1F *tigE_0p2aGate, *tigE_2p0aGate;
TH1F *tigE_0p2aDetSumGate, *tigE_2p0aDetSumGate, *tigE_horDetSumGate;
TH1F *tigE_2p0aDetSumGate_Downstream, *tigE_2p0aDetSumGate_Upstream, *tigE_2p0aDetSumGate_Corona;

TH1F *tigE_0p2aDetSumGate_Downstream, *tigE_0p2aDetSumGate_Upstream, *tigE_0p2aDetSumGate_Corona;

TH1F *tipChannel_HorizontalCut;
TH2F *tigE_tipRing;
TH2F *tigE_tipMult;



//PID gates
TCutG *alphaPosCut[NTIP], *protonPosCut[NTIP], *alphaPosCutSum, *protonPosCutSum, *horPosCutSum;

