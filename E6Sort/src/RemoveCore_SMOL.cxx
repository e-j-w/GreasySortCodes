//Sort code to remove crystals from SMOL files

#define RemoveCore_SMOL_cxx
#include "common.h"
#include "lin_eq_solver.h"
#include "evt_fmt.h"
#include "RemoveCore_SMOL.h"

using namespace std;

FILE *inp, *out;

void SortData(const char *sfile, const char *outfile, const int remCore[MAX_INPUT_CORE], const int numCoreVals){

    if((numCoreVals > MAX_INPUT_CORE)||(numCoreVals < 1)){
        cout << "ERROR: invalid number of input core values." << endl;
        exit(-1);
    }
    
    uint64_t sentries = 0U;
    uint64_t pileupCtrs[16];
    sorted_evt sortedEvt, sortedEvtOut;
    uint64_t hitBuildFlags = 0;
    uint64_t numSeparatedEvents = 0;

    cout << endl << "Removing cores from file: " << sfile << endl;
    cout << "Will write results to file: " << outfile << endl;

    //open the file to be sorted from
    inp = fopen(sfile, "r");
    if(inp == NULL){
        cout << "ERROR: couldn't open file " << sfile << endl;
        return;
    }
    fread(&sentries,sizeof(uint64_t),1,inp);
    uint64_t smolVersion = (uint64_t)(sentries >> 48);
    if(smolVersion > 0){
        fread(&pileupCtrs,sizeof(pileupCtrs),1,inp);
    }
    
    //setup the output file
    out = fopen(outfile, "wb");
    if(out == NULL){
        cout << "ERROR: couldn't open output file " << outfile << endl;
        return;
    }
    fwrite(&sentries,sizeof(uint64_t),1,out);
    if(smolVersion > 0){
        fwrite(&pileupCtrs,sizeof(pileupCtrs),1,out);
    }
    sentries &= 0xFFFFFFFFFFFF; // only first 48 bits specify number of events

    uint64_t actualSepEntries = 0;
    for(Long64_t jentry = 0; jentry < sentries; jentry++){

        //read event from input file
        if(readSMOLEvent(inp,&sortedEvt)==0){
            cout << "ERROR: bad event data in entry " << jentry << "." << endl;
            exit(-1);
        }
        memset(&sortedEvtOut,0,sizeof(sortedEvtOut));
        memcpy(&sortedEvtOut.header,&sortedEvt.header,sizeof(evt_header));
        sortedEvtOut.header.numNoABHits = 0;

        //remove cores
        for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
            uint8_t corePos = sortedEvt.noABHit[noABHitInd].core & 63U;
            if(corePos < (NGRIFPOS*4)){
                uint8_t keep = 1;
                for(int i=0; i<numCoreVals; i++){
                    if(corePos == remCore[i]){
                        keep = 0;
                        break;
                    }
                }
                if(keep){
                    memcpy(&sortedEvtOut.noABHit[sortedEvtOut.header.numNoABHits],&sortedEvt.noABHit[noABHitInd],sizeof(hpge_hit));
                    sortedEvtOut.header.numNoABHits++;
                    if(sortedEvtOut.header.numNoABHits >= MAX_EVT_HIT){
                        break; //cannot fit any more hits
                    }
                }
            }
        }

        //write out data, if there's any left
        if(sortedEvtOut.header.numNoABHits > 0){
            fwrite(&sortedEvtOut.header,sizeof(evt_header),1,out);
            //write hits, without segment data (no segments for GRIFFIN)
            for(int i = 0; i<sortedEvtOut.header.numNoABHits;i++){
                fwrite(&sortedEvtOut.noABHit[i].timeOffsetNs,sizeof(float),1,out);
                fwrite(&sortedEvtOut.noABHit[i].energy,sizeof(float),1,out);
                fwrite(&sortedEvtOut.noABHit[i].tsDiff,sizeof(uint8_t),1,out);
                fwrite(&sortedEvtOut.noABHit[i].core,sizeof(uint8_t),1,out);
            }
            actualSepEntries++;
        }
        

        if(jentry % 9713 == 0)
            cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;

    }
    cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;

    //write the number of separated events to the beginning of the file
    fseek(out,0,SEEK_SET);
    actualSepEntries &= 0xFFFFFFFFFFFF; // only first 48 bits specify number of events
    uint64_t smolFormatVersion = 1;
    actualSepEntries |= (smolFormatVersion << 48);
    fwrite(&actualSepEntries,sizeof(uint64_t),1,out);
    printf("Wrote %lu separated events to: %s\n",actualSepEntries & 0xFFFFFFFFFFFF,outfile);
    fclose(out);

    fclose(inp);

    return;

}
int main(int argc, char **argv){

    char const *sfile;
    char const *soutfile;
    char outName[256];
    int remCore[MAX_INPUT_CORE];
    int numCoreVals = 0;

    if(argc <= 1){
        cout << "Arguments: RemoveCore_SMOL smol_file output_smolfile_suffix core1 core2 core3..." << endl;
        cout << "A code for removing crystals from SMOL trees (eg. to exclude bad data from HV trips)." << endl;
        cout << "  *smol_file* can be a single SMOL tree (extension .smol), or a list of SMOL trees (extension .list, one filepath per line)." << endl;
        cout << "  *core1*, *core2*, etc are up to " << MAX_INPUT_CORE << " crystal positions (0-indexed) which are to be removed from the output file." << endl;
        return 0;
    }else if((argc >= 4)&&(argc < (4+MAX_INPUT_CORE))){
        sfile = argv[1];
        soutfile = argv[2];
        for(int arg=3; arg<argc; arg++){
            remCore[arg-3] = atoi(argv[arg]);
        }
        numCoreVals = argc-3;
    }else{
        printf("Incorrect arguments\nArguments: RemoveCore_SMOL smol_file output_smolfile_suffix core1 core2 core3...\n");
        return 0;
    }

    printf("Starting RemoveCore_SMOL code\n");

    if(strcmp(soutfile,"")==0){
        cout << "ERROR: output suffix cannot be empty." << endl;
        return 0;
    }
    cout << "Output file suffix: " << soutfile << endl;

    char filePrefix[256];

    const char *dot = strrchr(sfile, '.'); //get the file extension
    if(dot==NULL){
        cout << "ERROR: couldn't get SMOL tree or list file name." << endl;
        return 0;
    }

    printf("Will remove data belonging to cores: ");
    for(int i=0;i<numCoreVals;i++){
        printf(" %i",remCore[i]);
    }
    printf("\n");

    if(strcmp(dot + 1, "smol") == 0){
        strncpy(filePrefix,sfile,256);
        const char *tok = strtok(filePrefix,"."); //get the filename without the extension
        if(tok!=NULL){
            snprintf(outName,255,"%s_%s.smol",basename(tok),soutfile);
            //printf("Will write to file: %s\n",outName);
            SortData(sfile, outName, remCore, numCoreVals);
        }else{
            cout << "ERROR: improperly formatted filename: " << sfile << endl;
            return 0;
        }
    }else if(strcmp(dot + 1, "list") == 0){
        printf("SMOL tree list: %s\n", sfile);
        
        FILE *listfile;
        char str[256];

        if((listfile=fopen(sfile,"r"))==NULL){
            cout << "ERROR: Cannot open the list file: " << sfile << endl;
            return 0;
        }else{
            while(!(feof(listfile))){//go until the end of file is reached
                if(fgets(str,256,listfile)!=NULL){ //get an entire line
                    str[strcspn(str, "\r\n")] = 0;//strips newline characters from the string
                    strncpy(filePrefix,str,256);
                    const char *tok = strtok(filePrefix,"."); //get the filename without the extension
                    if(tok!=NULL){
                        snprintf(outName,255,"%s_%s.smol",basename(tok),soutfile);
                        //printf("Will write to file: %s\n",outName);
                        SortData(str, outName, remCore, numCoreVals);
                    }else{
                        cout << "ERROR: improperly formatted filename: " << str << endl;
                        return 0;
                    }
                }
            }
        }
    }else{
        cout << "ERROR: improper file extension for *smol_file* argument (should be .smol or .list)." << endl;
        return 0;
    }

    return 0;
}
