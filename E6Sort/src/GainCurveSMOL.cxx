//use the Makefile!

#define GainCurveS_cxx
#include "common.cxx"
#include "GainCurveSMOL.h"

using namespace std;

sorted_evt sortedEvt;

void GainCurveS::SortData(const char *sfile, const uint64_t startNumSec)
{
    
    FILE *inp = fopen(sfile, "rb");
    if(inp==NULL){
        printf("ERROR: couldn't open file: %s\n",sfile);
        exit(-1);
    }
    
    uint64_t sentries = 0U;
    fread(&sentries,sizeof(uint64_t),1,inp);
    
    uint8_t footerVal;

    /*cout << "TIGRESS positions: " << endl;
    for(int det=1;det<17;det++){
        for(int cry=0;cry<4;cry++){
        TVector3 pos = tigress->GetPosition(det, cry, 0, 110., false);
        cout << "det: " << det << ", cry: " << cry << ", position: [ " << pos.x() << " " << pos.y() << " " << pos.z() << " ]" << endl;
        }
    }*/

    printf("Sorting data from file %s, offset by %lu seconds.\n",sfile,startNumSec);
    for(Long64_t jentry = 0; jentry < sentries; jentry++){

        //read event
        if(readSMOLEvent(inp,&sortedEvt)==0){
            cout << "ERROR: bad event data in entry " << jentry << "." << endl;
            exit(-1);
        }
        
        for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
            if(sortedEvt.noABHit[noABHitInd].core < (NTIGPOS*4)){
                if(sortedEvt.noABHit[noABHitInd].energy > MIN_HPGE_EAB){
                    Double_t tSec = ((sortedEvt.header.evtTimeNs + sortedEvt.noABHit[noABHitInd].timeOffsetNs)/(1.0E9)) + (double)startNumSec;
                    hpgeE_time[sortedEvt.noABHit[noABHitInd].core]->Fill(tSec/60.0, sortedEvt.noABHit[noABHitInd].energy);
                }
            }
        }

        if (jentry % 10000 == 0)
            cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
    } // analysis tree

    cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
    cout << endl << "Event sorting complete" << endl;
    
    fclose(inp);
}

int main(int argc, char **argv)
{

  GainCurveS *mysort = new GainCurveS();

  const char *lsfile;
  const char *outfile;
  printf("Starting GainCurveSMOL\n");

  if(argc == 1){
    cout << "Code sorts diagnostic histograms for online HPGe data" << endl;
    cout << "Arguments: GainCurveSMOL list_file output_file" << endl;
    cout << "The list file should contain two columns (space-delimited) with the SMOL tree filenames and the times (in seconds) at the start of each run." << endl;
    return 0;
  }else if(argc == 3){
    lsfile = argv[1];
    outfile = argv[2];
    printf("List file: %s\nOutput file: %s\n", lsfile, outfile);
  }else{
    printf("ERROR: too many arguments!\nArguments: DecayCurve SMOL smol_file output_file\n");
    return 0;
  }

  theApp=new TApplication("App", &argc, argv);

  const char *dot = strrchr(lsfile, '.'); //get the file extension
  if(dot==NULL){
    cout << "ERROR: couldn't get list file name." << endl;
    return 0;
  }

  if(strcmp(dot + 1, "list") == 0){
    
    FILE *listfile;
    char str[256];
    char fileName[128];
    char *tok;

    if((listfile=fopen(lsfile,"r"))==NULL){
        cout << "ERROR: Cannot open the list file: " << lsfile << endl;
        return 0;
    }else{
        mysort->Initialise();
        while(!(feof(listfile))){//go until the end of file is reached
            if(fgets(str,256,listfile)!=NULL){ //get an entire line
                str[strcspn(str, "\r\n")] = 0;//strips newline characters from the string
                tok = strtok(str,"\t ");
                if(tok!=NULL){
                    strncpy(fileName,tok,127);
                    tok = strtok(NULL,"");
                    if(tok!=NULL){
                        mysort->SortData(fileName,(uint64_t)strtol(tok,NULL,10));
                    }else{
                        printf("ERROR: couldn't get time for run in file: %s\n",fileName);
                        exit(-1);
                    }
                }
            }
        }
    }
  }else{
    cout << "ERROR: improper file extension for input list (should be .list)." << endl;
    return 0;
  }

  TFile *myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  TDirectory *gcdir = myfile->mkdir("GainCurve");
  gcdir->cd();
  gcList->Write();
  myfile->cd();

  myfile->Write();
  myfile->Close();

  return 0;
}
