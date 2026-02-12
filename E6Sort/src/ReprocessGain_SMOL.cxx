//Sort code to re-do gain matching in SMOL files

#define ReprocessGain_SMOL_cxx
#include "common.h"
#include "lin_eq_solver.h"
#include "evt_fmt.h"
#include "ReprocessGain_SMOL.h"

using namespace std;

FILE *inp, *out;
double actualEnergy[MAX_INPUT_E];
double corrCoeff[MAX_INPUT_E][NGRIFPOS*4][MAX_EVT_WINDOWS_PER_TREE];
double enWindowAvg[MAX_INPUT_E][NGRIFPOS*4][MAX_EVT_WINDOWS_PER_TREE];
uint64_t enWindowNumHits[MAX_INPUT_E][NGRIFPOS*4][MAX_EVT_WINDOWS_PER_TREE];
uint8_t actualEnDetermined;

//fitter
long double xpowsum[5];//sums of (x1)^0, (x1)^1, (x1)^2, etc. indexed  by power #
long double mxpowsum[3];//sums of m*(x1)^0, m*(x1)^1, m*(x1)^2, etc. indexed by power #
  

void SortData(const char *sfile, const char *efile, const char *outfile, const uint64_t evalWindowSize, const uint64_t evtWindowSize, const double en[MAX_INPUT_E], const int numEnVals){

    if((numEnVals > MAX_INPUT_E)||(numEnVals < 1)){
        cout << "ERROR: invalid number of input energy values." << endl;
        exit(-1);
    }
    
    uint64_t sentries = 0U;
    uint64_t pileupCtrs[16];
    sorted_evt sortedEvt;
    uint64_t hitBuildFlags = 0;
    uint64_t numSeparatedEvents = 0;
    memset(enWindowAvg,0,sizeof(enWindowAvg));
    memset(enWindowNumHits,0,sizeof(enWindowNumHits));
    memset(corrCoeff,0,sizeof(corrCoeff));

    //evaluate proper gain, if it hasn't already been done
    if(actualEnDetermined == 0){

        inp = fopen(efile, "r");
        if(inp == NULL){
            cout << "ERROR: couldn't open file " << efile << endl;
            return;
        }
        fread(&sentries,sizeof(uint64_t),1,inp);
        uint64_t smolVersion = (uint64_t)(sentries >> 48);
        if(smolVersion > 0){
            fread(&pileupCtrs,sizeof(pileupCtrs),1,inp);
        }
        sentries &= 0xFFFFFFFFFFFF; // only first 48 bits specify number of events

        cout << endl << "Evaluating actual gain correction energies (from file " << efile << ")..." << endl;
        for(Long64_t jentry = 0; jentry < sentries; jentry++){

            //read event from input file
            if(readSMOLEvent(inp,&sortedEvt)==0){
                cout << "ERROR: bad event data in entry " << jentry << "." << endl;
                exit(-1);
            }

            if(jentry < evalWindowSize){
                for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
                    for(uint8_t i=0;i<numEnVals;i++){
                        if(fabs((sortedEvt.noABHit[noABHitInd].energy - en[i])/sortedEvt.noABHit[noABHitInd].energy) < EN_WINDOW_FRAC_WIDTH){
                            enWindowAvg[i][0][0] += sortedEvt.noABHit[noABHitInd].energy;
                            enWindowNumHits[i][0][0]++;
                        }
                    }
                }
            }else{
                break;
            }

            if(jentry % 9713 == 0)
                cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / evalWindowSize << "% complete" << "\r" << flush;
            
        }
        //calculate actual gain correction energies
        for(uint8_t i=0;i<numEnVals;i++){
            actualEnergy[i] = enWindowAvg[i][0][0]/(1.0*enWindowNumHits[i][0][0]);
        }
        printf("Determined acutal gain correction energies: ");
        for(uint8_t i=0;i<(numEnVals-1);i++){
            printf("%0.3f, ",actualEnergy[i]);
        }
        printf("%0.3f keV\n",actualEnergy[numEnVals-1]);
        actualEnDetermined = 1;
        //setup for next step
        memset(enWindowAvg,0,sizeof(enWindowAvg));
        memset(enWindowNumHits,0,sizeof(enWindowNumHits));

        fclose(inp); //done with energy evaluation file

        //open the file to be sorted
        inp = fopen(sfile, "r");
        if(inp == NULL){
            cout << "ERROR: couldn't open file " << sfile << endl;
            return;
        }
        fread(&sentries,sizeof(uint64_t),1,inp);
        smolVersion = (uint64_t)(sentries >> 48);
        if(smolVersion > 0){
            fread(&pileupCtrs,sizeof(pileupCtrs),1,inp);
        }
        sentries &= 0xFFFFFFFFFFFF; // only first 48 bits specify number of events
    }else{
        //open the file to be sorted
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
        sentries &= 0xFFFFFFFFFFFF; // only first 48 bits specify number of events
    }

    cout << endl << "Computing gain corrections for file: " << sfile << endl;
    cout << "Will write results to file: " << outfile << endl;

    //evaluate gain in window
    cout << endl << "Determining correction coefficients from data..." << endl;
    int numWindows = 0;
    for(Long64_t jentry = 0; jentry < sentries; jentry++){

        //read event from input file
        if(readSMOLEvent(inp,&sortedEvt)==0){
            cout << "ERROR: bad event data in entry " << jentry << "." << endl;
            exit(-1);
        }

        Double_t tSec = ((sortedEvt.header.evtTimeNs)/(1.0E9));

        //discard bad timing events
        if((tSec <= 0.0)||(tSec > MAX_RUN_TIME_SEC)){
            /*if(sortedEvt.header.numNoABHits > 0){
                printf("ERROR: zero time for event %li with hit data:\n",jentry);
                for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
                    printf("Hit %i, core: %u, energy: %f\n",noABHitInd,sortedEvt.noABHit[noABHitInd].core & 63U,sortedEvt.noABHit[noABHitInd].energy);
                    printf("tSec: %f\n",tSec);
                    //exit(-1);
                }
            }*/
            continue;
        }
        //printf("tSec: %f\n");

        //figure out what event window we're in
        int windowNum = (int)(floor(jentry/(1.0*evtWindowSize)));
        if(windowNum >= MAX_EVT_WINDOWS_PER_TREE){
            printf("WARNING: entry %li has invalid window number %i",jentry,windowNum);
            windowNum = MAX_EVT_WINDOWS_PER_TREE-1;
        }else if(windowNum < 0){
            printf("WARNING: entry %li has invalid window number %i",jentry,windowNum);
            windowNum = 0;
        }
        if((windowNum + 1) > numWindows){
            numWindows = windowNum + 1;
        }

        for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
            if((sortedEvt.noABHit[noABHitInd].core & 63U) < (NGRIFPOS*4)){
                for(uint8_t i=0;i<numEnVals;i++){
                    if(fabs((sortedEvt.noABHit[noABHitInd].energy - actualEnergy[i])/sortedEvt.noABHit[noABHitInd].energy) < EN_WINDOW_FRAC_WIDTH){
                        enWindowAvg[i][sortedEvt.noABHit[noABHitInd].core & 63U][windowNum] += sortedEvt.noABHit[noABHitInd].energy;
                        enWindowNumHits[i][sortedEvt.noABHit[noABHitInd].core & 63U][windowNum]++;
                    }
                }
            }else{
                printf("Invalid core for event %li hit %i: %u\n",jentry,noABHitInd,sortedEvt.noABHit[noABHitInd].core & 63U);
            }
        }

        if(jentry % 9713 == 0)
            cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
    }
    cout << "Entry " << sentries << " of " << sentries << ", 100% complete, " << numWindows << " gain correction windows processed." << endl;
    printf("Fitting...\n");
    for(int window = 0; window < numWindows; window++){
        //fit gain correction coefficients and copy them
        for(uint8_t coreNum=0; coreNum<(NGRIFPOS*4); coreNum++){

            //generate sums needed for fit
            memset(xpowsum,0,sizeof(xpowsum));
            memset(mxpowsum,0,sizeof(mxpowsum));
            for(uint8_t i=0;i<numEnVals;i++){
                double dataVal = (enWindowAvg[i][coreNum][window])/(1.0*enWindowNumHits[i][coreNum][window]);
                for(int sumPow = 0; sumPow < 5; sumPow ++){
                    xpowsum[sumPow] += pow(dataVal,sumPow*1.0);
                    if(sumPow < 3){
                        mxpowsum[sumPow] += actualEnergy[i]*pow(dataVal,sumPow*1.0);
                    }
                }
            }

            //construct equations (n=1 specific case)
            int i,j;
            long double fitPar[MAX_DIM]; //array holding parameters (desribing parboloid) from chisq minimization
            int lines;//number of data points
            lin_eq_type linEq;
            linEq.dim=3; //2nd order polynomial fit
            
            linEq.matrix[0][0]=xpowsum[4];
            linEq.matrix[0][1]=xpowsum[3];
            linEq.matrix[0][2]=xpowsum[2];
            
            linEq.matrix[1][1]=xpowsum[2];
            linEq.matrix[1][2]=xpowsum[1];
            
            linEq.matrix[2][2]=xpowsum[0];//bottom right entry
            
            //mirror the matrix (top right half mirrored to bottom left half)
            for(i=1;i<3;i++){
                for(j=0;j<i;j++){
                    linEq.matrix[i][j]=linEq.matrix[j][i];
                }
            }
            
            linEq.vector[0]=mxpowsum[2];
            linEq.vector[1]=mxpowsum[1];
            linEq.vector[2]=mxpowsum[0];
                
            //solve system of equations and assign values
            if(!(solve_lin_eq(&linEq)==1)){
                printf("ERROR: Could not determine fit parameters (par1).\n");
                printf("Perhaps there are not enough data points to perform a fit?\n");
                printf("Otherwise you can also try adjusting the fit range using the UPPER_LIMITS and LOWER_LIMITS options.\n");
                exit(-1);
            }
            
            //save fit parameters
            for(i=0;i<3;i++){
                corrCoeff[2-i][coreNum][window]=linEq.solution[i];
            }

            //if fit fails, set default values
            if((corrCoeff[0][coreNum][window] != corrCoeff[0][coreNum][window])||(corrCoeff[1][coreNum][window] != corrCoeff[1][coreNum][window])||(corrCoeff[2][coreNum][window] != corrCoeff[2][coreNum][window])){
                /*printf("Failed fit (core %u, window %u).\nEnergy values:\n",coreNum,window);
                for(uint8_t i=0;i<numEnVals;i++){
                    printf("%f %lu %f\n",enWindowAvg[i][coreNum][window],enWindowNumHits[i][coreNum][window],(enWindowAvg[i][coreNum][window])/(1.0*enWindowNumHits[i][coreNum][window]));
                }
                getc(stdin);*/
                corrCoeff[0][coreNum][window] = 0.0;
                corrCoeff[1][coreNum][window] = 0.0;
                corrCoeff[2][coreNum][window] = 0.0;
                
            }

        }
    }

    printf("\nCorrection coefficients:\n");
    for(int i=0; i<numWindows; i++){
        for(int j=0; j<(NGRIFPOS*4); j++){
            if((corrCoeff[0][j][i] != 0.0)&&(corrCoeff[1][j][i] != 0.0)){
                printf("Window %i, crystal %i: %0.4f %0.4f %0.4f\n",i,j,corrCoeff[0][j][i],corrCoeff[1][j][i],corrCoeff[2][j][i]);
            }
        }
    }

    //close and re-open the file
    fclose(inp);
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

    cout << endl << "Writing out corrected gains to file: " << outfile << endl;
    uint64_t actualSepEntries = 0;
    for(Long64_t jentry = 0; jentry < sentries; jentry++){

        //read event from input file
        if(readSMOLEvent(inp,&sortedEvt)==0){
            cout << "ERROR: bad event data in entry " << jentry << "." << endl;
            exit(-1);
        }

        //discard bad timing events
        Double_t tSec = ((sortedEvt.header.evtTimeNs)/(1.0E9));
        if((tSec <= 0.0)||(tSec > MAX_RUN_TIME_SEC)){
            continue;
        }

        //figure out what event window we're in
        int windowNum = (int)(floor(jentry/(1.0*evtWindowSize)));
        if(windowNum >= numWindows){
            windowNum = numWindows-1;
        }else if(windowNum < 0){
            windowNum = 0;
        }

        //correct energies
        for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
            uint8_t corePos = sortedEvt.noABHit[noABHitInd].core & 63U;
            if(corePos < (NGRIFPOS*4)){
                float newE = (float)(corrCoeff[0][corePos][windowNum] + corrCoeff[1][corePos][windowNum]*sortedEvt.noABHit[noABHitInd].energy + corrCoeff[2][corePos][windowNum]*sortedEvt.noABHit[noABHitInd].energy*sortedEvt.noABHit[noABHitInd].energy);
                sortedEvt.noABHit[noABHitInd].energy = newE;
            }else{
                sortedEvt.noABHit[noABHitInd].energy = 0.0f;
            }
        }

        //write out data
        fwrite(&sortedEvt.header,sizeof(evt_header),1,out);
        //write hits, without segment data (no segments for GRIFFIN)
        for(int i = 0; i<sortedEvt.header.numNoABHits;i++){
            fwrite(&sortedEvt.noABHit[i].timeOffsetNs,sizeof(float),1,out);
            fwrite(&sortedEvt.noABHit[i].energy,sizeof(float),1,out);
            fwrite(&sortedEvt.noABHit[i].tsDiff,sizeof(uint8_t),1,out);
            fwrite(&sortedEvt.noABHit[i].core,sizeof(uint8_t),1,out);
        }
        actualSepEntries++;

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
    printf("Wrote %lu separated events to: %s",actualSepEntries & 0xFFFFFFFFFFFF,outfile);
    fclose(out);

    fclose(inp);

    return;

}
int main(int argc, char **argv){

    char const *sfile, *efile;
    char const *soutfile;
    char outName[256];
    double corrEnergy[MAX_INPUT_E];
    uint64_t evalWindow, evtWindow;
    int numEnVals = 0;

    if(argc <= 1){
        cout << "Arguments: ReprocessGain_SMOL smol_file eval_smol_file eval_event_window corr_event_window output_smolfile_suffix energy1 energy2 energy3..." << endl;
        cout << "A code for re-aligning gains in SMOL trees." << endl;
        cout << "  *smol_file* can be a single SMOL tree (extension .smol), or a list of SMOL trees (extension .list, one filepath per line)." << endl;
        cout << "  *eval_smol_file* is a single SMOL tree (extension .smol) which is used to evaluate the correct gain." << endl;
        cout << "  *eval_event_window* is the size of the window at the start of the data in *eval_smol_file* that is used to evaluate the correct gain, in events.  This should correspond to the longest time period where gain shift is not observed." << endl;
        cout << "  *corr_event_window* is the size of the window used to evaluate gain shift, in events." << endl;
        cout << "  *energy1*, *energy2*, etc are up to " << MAX_INPUT_E << " approximate energies (in keV) corresponding to peaks used to the fit. At least 3 energies must be specified (we are doing a quadratic fit)." << endl;
        return 0;
    }else if((argc >= 9)&&(argc < (6+MAX_INPUT_E))){
        sfile = argv[1];
        efile = argv[2];
        evalWindow = atoll(argv[3]);
        evtWindow = atoll(argv[4]);
        soutfile = argv[5];
        for(int arg=6; arg<argc; arg++){
            corrEnergy[arg-6] = atof(argv[arg]);
        }
        numEnVals = argc-6;
    }else{
        printf("Incorrect arguments\nArguments: ReprocessGain_SMOL smol_file eval_event_window corr_event_window output_smolfile_suffix energy1 energy2 energy3...\n");
        return 0;
    }

    printf("Starting ReprocessGain_SMOL code\n");
    if(strcmp(soutfile,"")==0){
        cout << "ERROR: output suffix cannot be empty." << endl;
        return 0;
    }

    if((corrEnergy[0] <= 0.0)||(corrEnergy[1] <= 0.0)||(corrEnergy[2] == 0.0)){
        cout << "ERROR: gain correction energies must all be greater than zero." << endl;
        return 0;
    }
    if((corrEnergy[0] == corrEnergy[1])||(corrEnergy[0] == corrEnergy[2])||(corrEnergy[2] == corrEnergy[1])){
        cout << "ERROR: gain correction energies must not be equal." << endl;
        return 0;
    }

    //swap energies
    if(corrEnergy[0] > corrEnergy[1]){
        double swap = corrEnergy[1];
        corrEnergy[1] = corrEnergy[0];
        corrEnergy[0] = swap;
    }
    if(corrEnergy[1] > corrEnergy[2]){
        double swap = corrEnergy[2];
        corrEnergy[2] = corrEnergy[1];
        corrEnergy[1] = swap;
    }
    if(corrEnergy[0] > corrEnergy[1]){
        double swap = corrEnergy[1];
        corrEnergy[1] = corrEnergy[0];
        corrEnergy[0] = swap;
    }

    cout << "Gain evaluation window: " << evalWindow << " events" << endl;
    cout << "Gain correction window: " << evtWindow << " events" << endl;
    printf("Gain correction energies: ");
    for(int i=0;i<(numEnVals-1);i++){
        printf("%0.3f, ",corrEnergy[i]);
    }
    printf("%0.3f keV\n",corrEnergy[numEnVals-1]);
    cout << "Output file suffix: " << soutfile << endl;

    memset(actualEnergy,0,sizeof(actualEnergy)); //zero out
    actualEnDetermined = 0;

    char filePrefix[256];

    const char *dot = strrchr(sfile, '.'); //get the file extension
    if(dot==NULL){
        cout << "ERROR: couldn't get SMOL tree or list file name." << endl;
        return 0;
    }

    const char *dote = strrchr(efile, '.'); //get the file extension
    if(dote==NULL){
        cout << "ERROR: couldn't get eval SMOL tree file name." << endl;
        return 0;
    }

    if(strcmp(dote + 1, "smol") == 0){
        if(strcmp(dot + 1, "smol") == 0){
            strncpy(filePrefix,sfile,256);
            const char *tok = strtok(filePrefix,"."); //get the filename without the extension
            if(tok!=NULL){
                snprintf(outName,255,"%s_%s.smol",basename(tok),soutfile);
                //printf("Will write to file: %s\n",outName);
                SortData(sfile, efile, outName, evalWindow, evtWindow, corrEnergy, numEnVals);
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
                            SortData(str, efile, outName, evalWindow, evtWindow, corrEnergy, numEnVals);
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
    }else{
        cout << "ERROR: improper file extension for *eval_smol_file* argument (should be .smol)." << endl;
        return 0;
    }
    

    

    return 0;
}
