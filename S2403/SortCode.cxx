//g++ SortCode.cxx -std=c++0x -I$GRSISYS/include -L$GRSISYS/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -L$GRSISYS/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT -o SortData

#include "SortCode.h"		 
void MakeTigress(TTigress& tigress) {

  for(int i=0;i<tigress.GetMultiplicity();i++) {
    
    TTigressHit* tig_hit = tigress.GetTigressHit(i);
    int num = tig_hit->GetArrayNumber();
    
    double en = tig_hit->GetEnergy();
    if(en < TIG_THRESH)
      continue;
    
    int seg_mult = tig_hit->GetNSegments();
    for(int j=0;j<seg_mult;j++) {
      
      TDetectorHit seg_hit = tig_hit->GetSegmentHit(j);
      double seg_en = seg_hit.GetEnergy();
      //double seg_ts = seg_hit.GetTimeStamp();
      int seg = seg_hit.GetSegment();
      int seg_num = num*8 + seg - 1;
      
      FillHist("Tigress","Segment_Summary",512,0,512,seg_num,2000,0,4000,seg_en);
      
    }

    if(BadSeg(tig_hit))
	    continue;	
    
    long tm = tig_hit->GetTime();
    bool sup = tig_hit->BGOFired();
    TVector3 pos = tig_hit->GetPosition();

    FillHist("Tigress","Tigress_Energy",8000,0,8000,en);
    FillHist("Tigress","Tigress_Summary",64,0,64,num,4000,0,4000,en);
    FillHist("Tigress","Tigress_Position",180,0,180,pos.Theta()*r2d,200,-200,200,pos.Phi()*r2d);
    FillHist("Tigress","Tigress_EvT",1000,0,4000,tm/1000000000.,2000,0,2000,en);

    FillHist("Tigress/Dets",Form("tigTime%02d",num),1000,0,4000,tm/1000000000.);
    FillHist("Tigress/Dets",Form("tigEnergy%02d",num),4000,0,4000,en);
    //FillHist("Tigress/Dets",Form("tigET%02d",num),1000,0,5000,tm/1000000000.,1000,0,4000,en);
    
    if(!sup) {
      FillHist("Tigress","Suppressed_Energy",8000,0,8000,en);
      FillHist("Tigress","Suppressed_Summary",64,0,64,num,4000,0,4000,en);
    }
    
    for(int j=i+1;j<tigress.GetMultiplicity();j++) {

      TTigressHit* tig_hit2 = tigress.GetTigressHit(j);
	  
      double en2 = tig_hit2->GetEnergy();
      if(en2 < TIG_THRESH)
	      continue;

      //int num2 = tig_hit2->GetArrayNumber();
      long tm2 = tig_hit2->GetTime();
      //TVector3 pos2 = tig_hit2->GetPosition();

      double tdiff = tm-tm2;
      FillHist("Tigress","Tigress_TimeEnergy",3000,-3000,3000,tdiff,1200,0,2000,en);

    } //end second tigress loop (j)
  } //end first tigress loop (i)

  for(int i=0;i<tigress.GetAddbackMultiplicity();i++) {

    TTigressHit* add_hit = tigress.GetAddbackHit(i);

    double add_en = add_hit->GetEnergy();
    if(add_en < TIG_THRESH)
      continue;

    int add_num = add_hit->GetArrayNumber();
    long add_tm = add_hit->GetTime();
    //TVector3 add_pos = add_hit->GetPosition();

    FillHist("Tigress","Addback_Energy",8000,0,8000,add_en);
    FillHist("Tigress","Addback_Summary",64,0,64,add_num,4000,0,4000,add_en);

    bool supAdd = add_hit->BGOFired();    
    if(!supAdd) {
      FillHist("Tigress","Suppressed_Addback_Energy,",8000,0,8000,add_en);
      FillHist("Tigress","Suppressed_Addback_Summary",64,0,64,add_num,8000,0,8000,add_en);
    }

  } //end tigress addback loop
  tigress.ResetAddback(); //see https://github.com/GRIFFINCollaboration/GRSISort/issues/1527
  
  return;
} //End MakeTigress()

void MakeS3(TS3& s3) {
  
  for(int i=0;i<s3.GetRingMultiplicity();i++) {
	
    TS3Hit* ring_hit = s3.GetRingHit(i);
    
    int ring_det = ring_hit->GetDetector();
    int ring = ring_hit->GetRing();
    int num = 56*(ring_det-1) + ring + 32;

    double ring_ch = ring_hit->GetCharge();
    double ring_en = ring_hit->GetEnergy();
    long ring_tm = ring_hit->GetTime();

    FillHist("S3","charge_summary_chan",112,0,112,num,5000,0,20000,ring_ch);
    FillHist("S3","s3_summary_chan",112,0,112,num,5000,0,6000,ring_en);
    FillHist("S3/Segs",Form("s3Time%03d",num),1000,0,4000,ring_tm/1000000000.);
    FillHist("S3","sumc",2e3,0,2e3,ring_hit->GetChannelNumber(), 1e3,0,1e3,ring_hit->GetCharge());	
  } //end s3 ring loop

  for(int i=0;i<s3.GetSectorMultiplicity();i++) {
	
    TS3Hit* sec_hit = s3.GetSectorHit(i);
    
    int sec_det = sec_hit->GetDetector();
    int sec = sec_hit->GetSector();
    int num = 56*(sec_det-1) + sec;

    double sec_ch = sec_hit->GetCharge();
    double sec_en = sec_hit->GetEnergy();
    long sec_tm = sec_hit->GetTime();

    FillHist("S3","charge_summary_chan",112,0,112,num,5000,0,20000,sec_ch);
    FillHist("S3","s3_summary_chan",112,0,112,num,5000,0,6000,sec_en);
    FillHist("S3/Segs",Form("s3Time%03d",num),1000,0,4000,sec_tm/1000000000.);
    FillHist("S3","sumc",2e3,0,2e3,sec_hit->GetChannelNumber(), 1e3,0,1e3,sec_hit->GetCharge());	
    
    for(int j=0;j<s3.GetRingMultiplicity();j++) {
	  
      TS3Hit* ring_hit = s3.GetRingHit(j);
      
      int ring_det = ring_hit->GetDetector();
      if(ring_det != sec_det)
	      continue;
      
      int ring = ring_hit->GetRing();
      double ring_en = ring_hit->GetEnergy();
      long ring_tm = ring_hit->GetTime();
      double tdiff = ring_tm - sec_tm;
      
      FillHist("S3",Form("ring_sec_TDiff_chan_Det%02d",sec_det),3000,-3000,3000,tdiff);

      FillHist("S3",Form("ring_sec_E_EDiff_chan_Det%02d",sec_det),
	       1024,0,512,sec_en,1024,-64,64,sec_en-ring_en);
      FillHist("S3",Form("ring_sec_E_TDiff_chan_Det%02d",sec_det),
	       1024,0,512,sec_en,3000,-3000,3000,tdiff);
      
      FillHist("S3",Form("sector_v_ring_chan_Det%02d",sec_det),24,0,24,ring,32,0,32,sec);
      FillHist("S3",Form("sectorEn_v_ringEn_chan_Det%02d",sec_det),
	       2500,0,500,ring_en,2500,0,500,sec_en);
      
      FillHist("S3",Form("Chan_S3_TDiff_Det%02d",sec_det),64,0,64,ring,3000,-3000,3000,tdiff);
      FillHist("S3",Form("Chan_S3_TDiff_Det%02d",sec_det),64,0,64,sec+32,3000,-3000,3000,tdiff);
      
    } //end ring loop
  } //end sector loop
      
  for(int i=0;i<s3.GetPixelMultiplicity();i++) {
	
    TS3Hit* s3_hit = s3.GetPixelHit(i);

    int det = s3_hit->GetDetector();
    int ring = s3_hit->GetRing();
    int sec = s3_hit->GetSector();

    int ring_num = 56*(det-1) + ring + 32;
    int sec_num = 56*(det-1) + sec;

    double s3_en = s3_hit->GetEnergy();
    long s3_tm = s3_hit->GetTime();
    TVector3 posS3 = GetS3Position(det,ring,sec);
//    posS3.SetX(posS3.X() + 1.0);	// calculated 4/09/23
//    posS3.SetY(posS3.Y() - 0.5);	// calculated 4/09/23

//    posS3.SetX(posS3.X() + 1.15);	//Calculated from run 59076
//    posS3.SetY(posS3.Y() - 1.25);	//Calculated from run 59076

    posS3.SetX(posS3.X() + 1.35);	//Calculated from run 59104
    posS3.SetY(posS3.Y() - 0.65);	//Calculated from run 59104

    //double theta = posS3.Theta()*r2d + gRandom->Uniform(-0.25,0.25);
    //double phi = posS3.Phi()*r2d + gRandom->Uniform(-11.25,11.25);

    double theta = posS3.Theta()*r2d;
    double phi = posS3.Phi()*r2d;

    FillHist("S3","s3_summary",112,0,112,sec_num,5000,0,500,s3_en);
    FillHist("S3","s3_summary",112,0,112,ring_num,5000,0,500,s3_en);

    FillHist("S3",Form("sector_v_ring_Det%02d",det),24,0,24,ring,32,0,32,sec);
    FillHist("S3",Form("s3_energy_v_ring_Det%02d",det),24,0,24,ring,5000,0,500,s3_en);

    TVector3 posS3Smear = GetS3Position(det,ring,sec,false,true);

//    posS3Smear.SetX(posS3Smear.X() + 1.0);	//Calculated 4/09/23
//    posS3Smear.SetY(posS3Smear.Y() - 0.5);	//Calculated 4/09/23

//    posS3Smear.SetX(posS3Smear.X() + 1.15);	//Calculated from run 59076
//    posS3Smear.SetY(posS3Smear.Y() - 1.25);	//Calculated from run 59076

    posS3Smear.SetX(posS3Smear.X() + 1.35);	//Calculated from run 59104
    posS3Smear.SetY(posS3Smear.Y() - 0.65);	//Calculated from run 59104

    FillHist("S3",Form("s3_XY_Det%02d",det),250,-40.,40.,posS3Smear.X(),250,-40.,40.,posS3Smear.Y());
    FillHist("S3",Form("s3_pos_Det%02d",det),360,-200,200,posS3Smear.Phi() * TMath::RadToDeg(),360,0.0,180.0, posS3Smear.Theta() * TMath::RadToDeg());
    FillHist("S3",Form("s3_rvp_Det%02d",det),
	     32,-31.*pi/32.,33.*pi/32.,posS3.Phi(),24,11.0,35.0,posS3.Perp());
    
    //Gated Spectra
    for(TCutG* gate : s3_gates) {

      const char* name = gate->GetName();
      const bool US = ((string)name).find("US") != string::npos;
      
      int gate_det = 2;
      if(US)
	gate_det = 1;

      if(det != gate_det)
	continue;

      if(!gate->IsInside(ring,s3_en))
	continue;

      FillHist("S3",Form("sector_v_ring_%s",name),24,0,24,ring,32,0,32,sec);
      FillHist("S3",Form("s3_energy_v_ring_%s",name),24,0,24,ring,5000,0,500,s3_en);

      FillHist("S3",Form("s3_XY_%s",name),
	       250,-40.,40.,posS3.X(),250,-40.,40.,posS3.Y());
      //FillHist("S3",Form("s3_pos_%s",name),360,-200,200,phi,360,0.0,180.0,theta);
      FillHist("S3",Form("s3_rvp_%s",name),
	       32,-31.*pi/32.,33.*pi/32.,posS3.Phi(),24,11.0,35.0,posS3.Perp());
      
    } //End loop over S3 gates
  } //end s3 pixel loop
  
  return;
} //End MakeS3()

void MakeS3Tigress(TS3& s3, TTigress& tigress) {

  Ep = 335.;

  for(int i=0;i<s3.GetSectorMultiplicity();i++) {
	
    TS3Hit* s3_hit = s3.GetSectorHit(i);
    long s3_tm = s3_hit->GetTime();

    int sec_det = s3_hit->GetDetector();
    int sec = s3_hit->GetSector();
    int s3_num = 56*(sec_det-1) + sec;
    
    for(int j=0;j<tigress.GetMultiplicity();j++) {
      
      TTigressHit* tig_hit = tigress.GetTigressHit(j);
      int tig_num = tig_hit->GetArrayNumber();
      long tig_tm = tig_hit->GetTime();
      double tdiff = tig_tm - s3_tm;

      double tig_en = tig_hit->GetEnergy();
      
      FillHist("S3_Tigress","Tigress_S3sec_TE_chan",3000,-3000,3000,tdiff,4000,0,2000,tig_en);
      FillHist("S3_Tigress","TDiff_S3summary_chan",112,0,112,s3_num,3000,-3000,3000,tdiff);
    }

  }

  for(int i=0;i<s3.GetRingMultiplicity();i++) {
	
    TS3Hit* s3_hit = s3.GetRingHit(i);
    long s3_tm = s3_hit->GetTime();

    int ring_det = s3_hit->GetDetector();
    int ring = s3_hit->GetRing();
    int s3_num = 56*(ring_det-1) + ring + 32;

    for(int j=0;j<tigress.GetMultiplicity();j++) {
      
      TTigressHit* tig_hit = tigress.GetTigressHit(j);
      int tig_num = tig_hit->GetArrayNumber();
      long tig_tm = tig_hit->GetTime();
      double tdiff = tig_tm - s3_tm;

      double tig_en = tig_hit->GetEnergy();
      
      FillHist("S3_Tigress","Tigress_S3ring_TE_chan",3000,-3000,3000,tdiff,4000,0,2000,tig_en);
      FillHist("S3_Tigress","TDiff_S3summary_chan",112,0,112,s3_num,3000,-3000,3000,tdiff);
    }

  }
  
  for(int i=0;i<s3.GetPixelMultiplicity();i++) {
	
    TS3Hit* s3_hit = s3.GetPixelHit(i);

    int ring = s3_hit->GetRing();
    double s3_en = s3_hit->GetEnergy();
    
    int s3_det = s3_hit->GetDetector();
    int sec = s3_hit->GetSector();
	
    long s3_tm = s3_hit->GetTime();

    //TVector3 s3_pos = s3_hit->GetPosition();
    TVector3 s3_pos = GetS3Position(s3_det,ring,sec);
//    s3_pos.SetX(s3_pos.X() + 1.0);	
//    s3_pos.SetY(s3_pos.Y() - 0.5);	

//    s3_pos.SetX(s3_pos.X() + 1.15);	//Calculated from run 59076
//    s3_pos.SetY(s3_pos.Y() - 1.25);	//Calculated from run 59076

    s3_pos.SetX(s3_pos.X() + 1.35);	//Calculated from run 59104
    s3_pos.SetY(s3_pos.Y() - 0.65);	//Calculated from run 59104

    //Ungated spectra
    for(int j=0;j<tigress.GetMultiplicity();j++) {
      
      TTigressHit* tig_hit = tigress.GetTigressHit(j);
      
      double tig_en = tig_hit->GetEnergy();
      if(tig_en < TIG_THRESH)
	continue;

      if(BadSeg(tig_hit))
	continue;
      
      //int tig_num = tig_hit->GetArrayNumber();
      //double tig_dop = tig_hit->GetDoppler(beta);
      //TVector3 tig_pos = tig_hit->GetPosition();

      long tig_tm = tig_hit->GetTime();
      double tdiff = tig_tm - s3_tm;      

      FillHist("S3_Tigress","Tigress_S3_TDiff",3000,-3000,3000,tdiff);
      FillHist("S3_Tigress","Tigress_S3_TE",3000,-3000,3000,tdiff,2000,0,4000,tig_en);

      if(!Gate1D(tdiff,tigS3T[0],tigS3T[1]))
	continue;      

      FillHist("S3_Tigress",Form("Tigress_Energy_sDet%02d",s3_det),8000,0,8000,tig_en);
      
    } //end tigress loop
    
    //Gated Spectra
    for(TCutG* gate : s3_gates) {

      const char* name = gate->GetName();
      const bool US = ((string)name).find("US") != string::npos;
      
      int gate_det = 2;
      if(US)
	gate_det = 1;

      if(s3_det != gate_det)
	continue;
      
      if(!gate->IsInside(ring,s3_en))
	continue;
      
      //Variables for Doppler correction
      //These change based on what you're gating on
      //Calculate them based on the gate
      double beta, gam, recon_beta, recon_gam; //Velocities
      TVector3 rPos; //Reconstructed position
      
      const bool proj = ((string)name).find("Rb") != string::npos;
      SetupKinematics(proj,s3_pos,beta,gam,recon_beta,recon_gam,rPos);
      
      for(int j=0;j<tigress.GetMultiplicity();j++) {
	
	TTigressHit* tig_hit = tigress.GetTigressHit(j);

	double tig_en = tig_hit->GetEnergy();
	if(tig_en < TIG_THRESH)
	  continue;

	long tig_tm = tig_hit->GetTime();
	double tdiff = tig_tm - s3_tm;
	
	if(!Gate1D(tdiff,tigS3T[0],tigS3T[1]))
	  continue;

	if(BadSeg(tig_hit))
	  continue;

	int tig_num = tig_hit->GetArrayNumber();
	
	TVector3 tig_pos = tig_hit->GetPosition();
	double theta = s3_pos.Angle(tig_pos);
	double dop_en = gam*(1.0 - beta*TMath::Cos(theta))*tig_en;

	TVector3 reacPlane = s3_pos.Cross(TVector3(0.0,0.0,1.0));
	TVector3 detPlane = tig_pos.Cross(TVector3(0.0,0.0,1.0));

	double reac_phi = reacPlane.Phi();
	double det_phi = detPlane.Phi();

	double planeAng = reac_phi - det_phi;
	if(planeAng < 0)
	  planeAng += TMath::TwoPi();

	double recon_theta = rPos.Angle(tig_pos);
	double rec_en = recon_gam*(1.0 - recon_beta*TMath::Cos(recon_theta))*tig_en;

	TVector3 reconPlane = rPos.Cross(TVector3(0.0,0.0,1.0));
	double recon_phi = reconPlane.Phi();

	double reconAng = recon_phi - det_phi;
	if(reconAng < 0)
	  reconAng += TMath::TwoPi();

	FillHist(Form("S3_Tigress_%s",name),Form("Tigress_Energy_%s",name),
		 4000,0,4000,tig_en);
	FillHist(Form("S3_Tigress_%s",name),Form("Doppler_Energy_%s",name),
		 4000,0,4000,dop_en);
	FillHist(Form("S3_Tigress_%s",name),Form("Recon_Energy_%s",name),
		 4000,0,4000,rec_en);

	if(tigress.GetMultiplicity() == 1){
	  FillHist(Form("S3_Tigress_%s",name),Form("Tigress_Energy_Mult1_%s",name),
	  	   4000,0,4000,tig_en);
	  FillHist(Form("S3_Tigress_%s",name),Form("Doppler_Energy_Mult1_%s",name),
		   4000,0,4000,dop_en);
	  FillHist(Form("S3_Tigress_%s",name),Form("Recon_Energy_Mult1_%s",name),
		   4000,0,4000,rec_en);
	}

	if(tigress.GetMultiplicity() >= 2) {
	  FillHist(Form("S3_Tigress_%s",name),Form("Tigress_Energy_Mult2_%s",name),
	  	   4000,0,4000,tig_en);
	  FillHist(Form("S3_Tigress_%s",name),Form("Doppler_Energy_Mult2_%s",name),
		   4000,0,4000,dop_en);
	  FillHist(Form("S3_Tigress_%s",name),Form("Recon_Energy_Mult2_%s",name),
		   4000,0,4000,rec_en);
	}
	
	FillHist(Form("S3_Tigress_%s",name),Form("Tigress_Summay_%s",name),
		 64,0,64,tig_num,4000,0,4000,tig_en);
	FillHist(Form("S3_Tigress_%s",name),Form("Doppler_Summary_%s",name),
		 64,0,64,tig_num,4000,0,4000,dop_en);
	FillHist(Form("S3_Tigress_%s",name),Form("Recon_Summary_%s",name),
		 64,0,64,tig_num,4000,0,4000,rec_en);

	/*	
	if(proj) {
	  for(int i=0;i<360;i++) {

	    double phi_off = i + 0.5;
	    TVector3 tmp_pos(s3_pos);
	    tmp_pos.SetPhi(tmp_pos.Phi() + phi_off*d2r);

	    double tmp_theta = tmp_pos.Angle(tig_pos);
	    double tmp_dop = gam*(1.0 - beta*TMath::Cos(tmp_theta))*tig_en;

	    FillHist(Form("S3_Tigress_%s",name),Form("Phi_Scan_%s",name),
		     360,0,360,phi_off,4000,0,4000,tmp_dop);
	    
	  } 
	}
	*/	

	FillHist(Form("GammaCorrelations_%s",name),Form("Theta_Correlation_%s",name),
		 2000,0,2000,tig_en,90,0,180,theta*r2d);
	FillHist(Form("GammaCorrelations_%s",name),Form("Theta_Correction_%s",name),
		 2000,0,2000,dop_en,90,0,180,theta*r2d);
	


	/*for(int i=0;i<21;i++){
		Ep	= 260 + 4. * i;
		SetupKinematics(proj,s3_pos,beta,gam,recon_beta,recon_gam,rPos);

		double dop_en2 = gam*(1.0 - beta*TMath::Cos(theta))*tig_en;
		double rec_en_mod = recon_gam*(1.0 - recon_beta*TMath::Cos(recon_theta))*tig_en;

		FillHist(Form("GammaCorrelations_%s/IterateCorrection/",name),Form("Projection_Theta_Correlation_DopEn_%s_%f",name,Ep),
			 2000,0,2000,dop_en2);
		FillHist(Form("GammaCorrelations_%s/IterateCorrection/",name),Form("Projection_Theta_Correction_RecEn_%s_%f",name,Ep),
			 2000,0,2000,rec_en_mod);
	
		FillHist(Form("GammaCorrelations_%s/IterateCorrection/",name),Form("Theta_Correlation_DopEn_%s_%f",name,Ep),
			 2000,0,2000,dop_en2,90,0,180,theta*r2d);
		FillHist(Form("GammaCorrelations_%s/IterateCorrection/",name),Form("Theta_Correction_RecEn_%s_%f",name,Ep),
			 2000,0,2000,rec_en_mod,90,0,180,recon_theta*r2d);
		
		FillHist(Form("GammaCorrelations_%s/IterateCorrection/",name),Form("S3Theta_Correlation_DopEn_%s_%f",name,Ep),
			 2000,0,2000,dop_en2,90,0,180,s3_pos.Theta()*r2d);
		FillHist(Form("GammaCorrelations_%s/IterateCorrection/",name),Form("S3RecTheta_Correction_RecEn_%s_%f",name,Ep),
			 2000,0,2000,rec_en_mod,90,0,180,rPos.Theta()*r2d);
	
		FillHist(Form("GammaCorrelations_%s/IterateCorrection/",name),Form("Phi_Correlation_DopEn_%s_%f",name,Ep),
			 2000,0,2000,dop_en2,90,0,360,planeAng*r2d);
		FillHist(Form("GammaCorrelations_%s/IterateCorrection/",name),Form("Phi_Correction_RecEn_%s_%f",name,Ep),
			 2000,0,2000,rec_en_mod,90,0,360,reconAng*r2d);
						
	
	}*/
	
	FillHist(Form("GammaCorrelations_%s",name),Form("Phi_Correlation_%s",name),
		 2000,0,2000,tig_en,90,0,360,planeAng*r2d);
	FillHist(Form("GammaCorrelations_%s",name),Form("Phi_Correction_%s",name),
		 2000,0,2000,dop_en,90,0,360,planeAng*r2d);

	FillHist(Form("GammaCorrelations_%s",name),Form("recTheta_Correlation_%s",name),
		 2000,0,2000,tig_en,90,0,180,recon_theta*r2d);
	FillHist(Form("GammaCorrelations_%s",name),Form("recTheta_Correction_%s",name),
		 2000,0,2000,rec_en,90,0,180,recon_theta*r2d);

	FillHist(Form("GammaCorrelations_%s",name),Form("recPhi_Correlation_%s",name),
		 2000,0,2000,tig_en,90,0,360,reconAng*r2d);
	FillHist(Form("GammaCorrelations_%s",name),Form("recPhi_Correction_%s",name),
		 2000,0,2000,rec_en,90,0,360,reconAng*r2d);

	FillHist(Form("Angle_Ranges_%s",name),Form("tigEn_v_Ring_%s",name),
		 24,0,24,ring,4000,0,4000,tig_en);
	FillHist(Form("Angle_Ranges_%s",name),Form("dopEn_v_Ring_%s",name),
		 24,0,24,ring,4000,0,4000,dop_en);
	FillHist(Form("Angle_Ranges_%s",name),Form("recEn_v_Ring_%s",name),
		 24,0,24,ring,4000,0,4000,rec_en);

	if(ring < 6) {
	  FillHist(Form("Angle_Ranges_%s",name),Form("tigEn_%s_Y1",name),
		   4000,0,4000,tig_en);
	  FillHist(Form("Angle_Ranges_%s",name),Form("dopEn_%s_Y1",name),
		   4000,0,4000,dop_en);
	  FillHist(Form("Angle_Ranges_%s",name),Form("recEn_%s_Y1",name),
		   4000,0,4000,rec_en);
	}
	else if(ring < 12) {
	  FillHist(Form("Angle_Ranges_%s",name),Form("tigEn_%s_Y2",name),
		   4000,0,4000,tig_en);
	  FillHist(Form("Angle_Ranges_%s",name),Form("dopEn_%s_Y2",name),
		   4000,0,4000,dop_en);
	  FillHist(Form("Angle_Ranges_%s",name),Form("recEn_%s_Y2",name),
		   4000,0,4000,rec_en);
	}
	else if(ring < 18) {
	  FillHist(Form("Angle_Ranges_%s",name),Form("tigEn_%s_Y3",name),
		   4000,0,4000,tig_en);
	  FillHist(Form("Angle_Ranges_%s",name),Form("dopEn_%s_Y3",name),
		   4000,0,4000,dop_en);
	  FillHist(Form("Angle_Ranges_%s",name),Form("recEn_%s_Y3",name),
		   4000,0,4000,rec_en);
	}
	else {
	  FillHist(Form("Angle_Ranges_%s",name),Form("tigEn_%s_Y4",name),
		   4000,0,4000,tig_en);
	  FillHist(Form("Angle_Ranges_%s",name),Form("dopEn_%s_Y4",name),
		   4000,0,4000,dop_en);
	  FillHist(Form("Angle_Ranges_%s",name),Form("recEn_%s_Y4",name),
		   4000,0,4000,rec_en);
	}
	
	//FillHist(Form("S3_Tigress_%s",name),Form("_%s",name),);
	//FillHist(Form("S3_Tigress_%s",name),Form("_%s",name),);
	
      } //End tigress loop 

      for(int j=0;j<tigress.GetAddbackMultiplicity();j++) {
	
	TTigressHit* add_hit = tigress.GetAddbackHit(j);
	
	bool supAdd = add_hit->BGOFired();  
	if(supAdd)
	  continue;
	
	double add_en = add_hit->GetEnergy();
	if(add_en < TIG_THRESH)
	  continue;

	long tig_tm = add_hit->GetTime();
	double tdiff = tig_tm - s3_tm;
	
	if(!Gate1D(tdiff,tigS3T[0],tigS3T[1]))
	  continue;
	
	//int add_num = add_hit->GetArrayNumber();
	//long add_tm = add_hit->GetTime();

	
	TVector3 add_pos = add_hit->GetPosition();
	double add_theta = s3_pos.Angle(add_pos);

	double add_dop_en = gam*(1.0 - beta*TMath::Cos(add_theta))*add_en;

	double add_recon_theta = rPos.Angle(add_pos);
	double add_rec_en = recon_gam*(1.0 - recon_beta*TMath::Cos(add_recon_theta))*add_en;
	  
	FillHist(Form("S3_Tigress_%s",name),Form("Tigress_Addback_Energy_%s",name),
		 4000,0,4000,add_en);
	FillHist(Form("S3_Tigress_%s",name),Form("Doppler_Addback_Energy_%s",name),
		 4000,0,4000,add_dop_en);
	FillHist(Form("S3_Tigress_%s",name),Form("Recon_Addback_Energy_%s",name),
		 4000,0,4000,add_rec_en);

	if(tigress.GetAddbackMultiplicity() == 1){
	  FillHist(Form("S3_Tigress_%s",name),Form("Tigress_Addback_Energy_Mult1_%s",name),
	  	   4000,0,4000,add_en);
	  FillHist(Form("S3_Tigress_%s",name),Form("Doppler_Addback_Energy_Mult1_%s",name),
		   4000,0,4000,add_dop_en);
	  FillHist(Form("S3_Tigress_%s",name),Form("Recon_Addback_Energy_Mult1_%s",name),
		   4000,0,4000,add_rec_en);
	}

	if(tigress.GetAddbackMultiplicity() >= 2) {
	  FillHist(Form("S3_Tigress_%s",name),Form("Tigress_Addback_Energy_Mult2_%s",name),
	  	   4000,0,4000,add_en);
	  FillHist(Form("S3_Tigress_%s",name),Form("Doppler_Addback_Energy_Mult2_%s",name),
		   4000,0,4000,add_dop_en);
	  FillHist(Form("S3_Tigress_%s",name),Form("Recon_Addback_Energy_Mult2_%s",name),
		   4000,0,4000,add_rec_en);
	}

	for(int k=j+1;k<tigress.GetAddbackMultiplicity();k++) {

		TTigressHit* add_hit2 = tigress.GetAddbackHit(k);

		double add_en2 = add_hit2->GetEnergy();
		if(add_en2 < TIG_THRESH)
		  continue;

		long tig_tm2 = add_hit2->GetTime();
		double tdiff2 = tig_tm2 - s3_tm;
	
		if(!Gate1D(tdiff2,tigS3T[0],tigS3T[1]))
	  		continue;
	
		//int add_num = add_hit->GetArrayNumber();
		//long add_tm = add_hit->GetTime();

	
		TVector3 add_pos2 = add_hit2->GetPosition();
		double add_theta = s3_pos.Angle(add_pos2);

		double add_dop_en2 = gam*(1.0 - beta*TMath::Cos(add_theta))*add_en2;

		double add_recon_theta2 = rPos.Angle(add_pos2);
		double add_rec_en2 = recon_gam*(1.0 - recon_beta*TMath::Cos(add_recon_theta2))*add_en2;
	  
		FillHist(Form("S3_Tigress_%s",name),Form("Tigress_AddbackAddback_Energy_%s",name),
		 	1000,0,1000,add_en,1000,0,1000,add_en2);
		FillHist(Form("S3_Tigress_%s",name),Form("Tigress_AddbackAddback_Energy_%s",name),
		 	1000,0,1000,add_en2,1000,0,1000,add_en);
		FillHist(Form("S3_Tigress_%s",name),Form("Doppler_AddbackAddback_Energy_%s",name),
		 	1000,0,1000,add_dop_en,1000,0,1000,add_dop_en2);
		FillHist(Form("S3_Tigress_%s",name),Form("Doppler_AddbackAddback_Energy_%s",name),
		 	1000,0,1000,add_dop_en2,1000,0,1000,add_dop_en);
		FillHist(Form("S3_Tigress_%s",name),Form("DopplerTigress_AddbackAddback_Energy_%s",name),
		 	1000,0,1000,add_dop_en,1000,0,1000,add_en2);
		FillHist(Form("S3_Tigress_%s",name),Form("DopplerTigress_AddbackAddback_Energy_%s",name),
		 	1000,0,1000,add_dop_en2,1000,0,1000,add_en);
		FillHist(Form("S3_Tigress_%s",name),Form("Recon_AddbackAddback_Energy_%s",name),
		 	1000,0,1000,add_rec_en,1000,0,1000,add_rec_en2);
		FillHist(Form("S3_Tigress_%s",name),Form("Recon_AddbackAddback_Energy_%s",name),
		 	1000,0,1000,add_rec_en2,1000,0,1000,add_rec_en);

	} 
	  
      } //end tigress addback loop

    } //End loop over S3 gates
  } //End S3 loop
  
  return;
}

void SortData(char const *afile, char const *calfile, char const *outfile) {
  
  cout << "Reading calibration file: " << calfile << endl;
  if(TChannel::ReadCalFile(calfile) < 1) {
    cout << "No channels found in calibration file " << calfile << "!" << endl;
    return;
  }
  
  cout << "Opening analysis tree file: " << afile << endl;
  TFile *analysisfile = new TFile(afile, "READ");
  if(!analysisfile->IsOpen()) {
    cout << " Failed to open file" << endl;
    return;
  }
  
  TChain* AnalysisTree = (TChain *)analysisfile->Get("AnalysisTree");
  if(!AnalysisTree) {
    cout << "Could not find AnalysisTree in file " << afile << "!" << endl;
    return;
  }
  long num_entries = AnalysisTree->GetEntries();

  TTigress* tigress = NULL;
  if(AnalysisTree->FindBranch("TTigress"))
    AnalysisTree->SetBranchAddress("TTigress",&tigress);
  else
    cout << "Branch 'TTigress' not found! TTigress variable is NULL pointer" << endl;
  
  TS3* s3 = NULL;
  if(AnalysisTree->FindBranch("TS3")) {
    AnalysisTree->SetBranchAddress("TS3",&s3);
    s3->SetFrontBackTime(160);
    //cout << "Setting loose front-back time condition (testing)." << endl;
    //s3->SetFrontBackTime(500000000);
    s3->SetFrontBackEnergy(0.9);
    //cout << "Setting loose front-back energy condition (testing)." << endl;
    //s3->SetFrontBackEnergy(0.0001);
    s3->PreferenceSector(false);
    s3->SetMultiHit(false);
  }
  else
    cout << "Branch 'TS3' not found! TS3 variable is NULL pointer" << endl;


  //TRF* RF = NULL;
  //if(AnalysisTree->FindBranch("TRF")) {
    //AnalysisTree->SetBranchAddress("TRF",&RF);
  //}
  //else
    //cout << "Branch 'TRF' not found! TRF variable is NULL pointer" << endl;
  
  hists = new TList();
  //LoadGates(s3_gate_names);
  
  srimB.ReadEnergyLossFile(beam_srim_name.c_str());
  srimT.ReadEnergyLossFile(targ_srim_name.c_str());

  /*for(int i=0;i<2;i++){
    for(int j=0;j<10;j++){
      for(int k=0;k<11;k++){
        hIter[i][j][k] = new TH2F(Form("Iter_%i_%i_%i",i,j,k),Form("Iter_%i_%i_%i",i,j,k),1000,0,1000,90,0,180);
      }
    }
  }*/

  cout << "\nSorting " << num_entries << " analysis events..." << endl;
  for(long jentry=0;jentry<num_entries;jentry++) {
    AnalysisTree->GetEntry(jentry);
    
    if(tigress)
      MakeTigress(*tigress);

    if(s3)
      MakeS3(*s3);

    if(s3 && tigress)
      MakeS3Tigress(*s3,*tigress);
    
    if(!(jentry%10000))
      cout << " Event " << jentry << " (" << 100 * jentry / num_entries << "%)" << "\r"
		<< flush;
    
  } //end analysis tree loop
  
  cout << "\r Event " << num_entries << " (100%)" << "\nWriting histograms to "
	    << outfile << endl;

  WriteHists(outfile);
  
  cout << "Done!" << endl;
  
  return;
}

int main(int argc, char **argv) {

  cout << "Starting sortcode" << endl;
  if(argc == 1) {
    cout << "Please provide arguments" << endl;
    cout << "USAGE: SortData analysis_tree cal_file (optional) out_file (optional)" << endl;
    return 0;
  }

  string grsi_path = getenv("GRSISYS");
  if (grsi_path.length() > 0) {
    grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);
  TParserLibrary::Get()->Load();
  
  char const *afile = argv[1]; //Input analysis tree
  char const *calfile = "CalibrationFile.cal"; //Input calibration file
  char const *outfile = "Histograms.root"; //Output histogram file
 
  if(argc > 2)
    calfile = argv[2];
  if(argc > 3)
    outfile = argv[3];
  if(argc > 4) {
    cout << "Too many arguments" << endl;
    cout << "USAGE: SortData analysis_tree cal_file (optional) out_file (optional)" << endl;
    return 0;
  }

  cout << "Analysis file: " << afile << "\nCalibration file: " << calfile << "\nOutput file: "
	    << outfile << endl;

  SortData(afile,calfile,outfile);

  return 0;
}
