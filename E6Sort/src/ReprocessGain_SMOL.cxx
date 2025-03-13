//Sort code to re-do gain matching in SMOL files

#define ReprocessGain_SMOL_cxx
#include "common.h"
#include "lin_eq_solver.h"
#include "evt_fmt.h"
#include "ReprocessGain_SMOL.h"

using namespace std;

FILE *inp, *out;
double actualEnergy[MAX_INPUT_E];
double corrCoeff[MAX_INPUT_E][NTIGPOS*4][MAX_TIME_WINDOWS_PER_TREE];
double enWindowAvg[MAX_INPUT_E][NTIGPOS*4][MAX_TIME_WINDOWS_PER_TREE];
uint64_t enWindowNumHits[MAX_INPUT_E][NTIGPOS*4][MAX_TIME_WINDOWS_PER_TREE];
uint8_t actualEnDetermined;

//fitter
long double xpowsum[5];//sums of (x1)^0, (x1)^1, (x1)^2, etc. indexed  by power #
long double mxpowsum[3];//sums of m*(x1)^0, m*(x1)^1, m*(x1)^2, etc. indexed by power #
  

void ReprocessGain_SMOL::SortData(const char *sfile, const char *outfile, const double evalWindowSize, const double tWindowSize, const double en[MAX_INPUT_E], const int numEnVals){

    if((numEnVals > MAX_INPUT_E)||(numEnVals < 1)){
        cout << "ERROR: invalid number of input energy values." << endl;
        exit(-1);
    }

    inp = fopen(sfile, "r");
    if(inp == NULL){
        cout << "ERROR: couldn't open file " << sfile << endl;
        return;
    }
    
    uint64_t sentries = 0U;
    sorted_evt sortedEvt;
    uint8_t footerVal = 227U;
    uint64_t hitBuildFlags = 0;
    uint64_t numSeparatedEvents = 0;
    memset(enWindowAvg,0,sizeof(enWindowAvg));
    memset(enWindowNumHits,0,sizeof(enWindowNumHits));
    memset(corrCoeff,0,sizeof(corrCoeff));

    fread(&sentries,sizeof(uint64_t),1,inp);

    cout << endl << "Computing gain corrections for file: " << sfile << endl;
    cout << "Will write results to file: " << outfile << endl;

    //evaluate proper gain
    if(actualEnDetermined == 0){
        cout << endl << "Evaluating actual gain correction energies..." << endl;
        for(Long64_t jentry = 0; jentry < sentries; jentry++){

            //read event from input file
            if(readSMOLEvent(inp,&sortedEvt)==0){
                cout << "ERROR: bad event data in entry " << jentry << "." << endl;
                exit(-1);
            }

            Double_t tSec = ((sortedEvt.header.evtTimeNs)/(1.0E9));
            if((tSec > 0.0)&&(tSec < evalWindowSize)){
                for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
                    for(uint8_t i=0;i<numEnVals;i++){
                        if(fabs(sortedEvt.noABHit[noABHitInd].energy - en[i]) < 5.0){
                            enWindowAvg[i][0][0] += sortedEvt.noABHit[noABHitInd].energy;
                            enWindowNumHits[i][0][0]++;
                        }
                    }
                }
            }

            if(jentry % 9713 == 0)
                cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
            
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

        //close and re-open the file
        fclose(inp);
        inp = fopen(sfile, "r");
        if(inp == NULL){
            cout << "ERROR: couldn't open file " << sfile << endl;
            return;
        }
        fread(&sentries,sizeof(uint64_t),1,inp);
    }

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

        //get rid of zero-time empty events
        if(tSec <= 0.0){
            if(sortedEvt.header.numNoABHits == 0){
                continue;
            }else{
                printf("ERROR: zero time for event with hit data:\n");
                for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
                    printf("Hit %i, core: %u, energy: %f\n",noABHitInd,sortedEvt.noABHit[noABHitInd].core,sortedEvt.noABHit[noABHitInd].energy);
                    exit(-1);
                }
            }
        }

        int windowNum = (int)(floor((tSec)/(1.0*tWindowSize)));
        if(windowNum >= MAX_TIME_WINDOWS_PER_TREE){
            printf("WARNING: entry %li has invalid window number %i",jentry,windowNum);
            windowNum = MAX_TIME_WINDOWS_PER_TREE-1;
        }else if(windowNum < 0){
            printf("WARNING: entry %li has invalid window number %i",jentry,windowNum);
            windowNum = 0;
        }
        if((windowNum + 1) > numWindows){
            numWindows = windowNum + 1;
        }

        for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
            if(sortedEvt.noABHit[noABHitInd].core < (NTIGPOS*4)){
                for(uint8_t i=0;i<numEnVals;i++){
                    if(fabs(sortedEvt.noABHit[noABHitInd].energy - actualEnergy[i]) < 5.0){
                        enWindowAvg[i][sortedEvt.noABHit[noABHitInd].core][windowNum] += sortedEvt.noABHit[noABHitInd].energy;
                        enWindowNumHits[i][sortedEvt.noABHit[noABHitInd].core][windowNum]++;
                    }
                }
            }
        }

        if(jentry % 9713 == 0)
            cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
    }
    cout << "Entry " << sentries << " of " << sentries << ", 100% complete, " << numWindows << " gain correction windows processed." << endl;
    cout << endl << "Fitting..." << endl;
    for(int window = 0; window < numWindows; window++){
        //fit gain correction coefficients and copy them
        for(uint8_t coreNum=0; coreNum<(NTIGPOS*4); coreNum++){

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
            if(corrCoeff[0][coreNum][window] != corrCoeff[0][coreNum][window]){
                corrCoeff[0][coreNum][window] = 0.0;
            }
            if(corrCoeff[1][coreNum][window] != corrCoeff[1][coreNum][window]){
                corrCoeff[1][coreNum][window] = 1.0;
            }
            if(corrCoeff[2][coreNum][window] != corrCoeff[2][coreNum][window]){
                corrCoeff[2][coreNum][window] = 0.0;
            }

        }
    }

    printf("\nCorrection coefficients:\n");
    for(int i=0; i<numWindows; i++){
        for(int j=0; j<(NTIGPOS*4); j++){
            printf("Window %i, crystal %i: %0.4f %0.4f %0.4f\n",i,j,corrCoeff[0][j][i],corrCoeff[1][j][i],corrCoeff[2][j][i]);
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
    
    //setup the output file
    out = fopen(outfile, "wb");
    if(out == NULL){
        cout << "ERROR: couldn't open output file " << outfile << endl;
        return;
    }
    fwrite(&sentries,sizeof(uint64_t),1,out);

    cout << endl << "Writing out corrected gains to file: " << outfile << endl;
    uint64_t actualSepEntries = 0;
    for(Long64_t jentry = 0; jentry < sentries; jentry++){

        //read event from input file
        if(readSMOLEvent(inp,&sortedEvt)==0){
            cout << "ERROR: bad event data in entry " << jentry << "." << endl;
            exit(-1);
        }

        //figure out what time window we're in
        Double_t tSec = ((sortedEvt.header.evtTimeNs)/(1.0E9));
        if(tSec <= 0.0){
            continue;
        }
        int windowNum = (int)(floor((tSec)/(1.0*tWindowSize)));
        if(windowNum >= numWindows){
            windowNum = numWindows-1;
        }else if(windowNum < 0){
            windowNum = 0;
        }

        //correct energies
        for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
            uint8_t corePos = sortedEvt.noABHit[noABHitInd].core;
            if(corePos < (NTIGPOS*4)){
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
            fwrite(&sortedEvt.noABHit[i].core,sizeof(uint8_t),1,out);
        }
        //write footer value
        fwrite(&footerVal,sizeof(uint8_t),1,out);
        actualSepEntries++;

        if(jentry % 9713 == 0)
            cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;

    }
    cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;

    //write the number of separated events to the beginning of the file
    fseek(out,0,SEEK_SET);
    fwrite(&actualSepEntries,sizeof(uint64_t),1,out);
    cout << "Wrote " << actualSepEntries << " separated events to: " << outfile << endl;
    fclose(out);

    fclose(inp);

    return;

}
int main(int argc, char **argv){

    ReprocessGain_SMOL *mysort = new ReprocessGain_SMOL();

    char const *sfile;
    char const *soutfile;
    char outName[64];
    double corrEnergy[MAX_INPUT_E];
    double evalWindow, tWindow;
    int numEnVals = 0;

    if(argc <= 1){
        cout << "Arguments: ReprocessGain_SMOL smol_file eval_time_window corr_time_window output_smolfile_suffix energy1 energy2 energy3..." << endl;
        cout << "  *smol_file* can be a single SMOL tree (extension .smole6), or a list of SMOL trees (extension .list, one filepath per line)." << endl;
        cout << "  *eval_time_window* is the size of the window at the start of the data used to evaluate the correct gain, in seconds.  This should be the longest time period at the start of the data where gain shift is not observed." << endl;
        cout << "  *corr_time_window* is the size of the window used to evaluate gain shift, in seconds." << endl;
        cout << "  *energy1*, *energy2*, etc are up to " << MAX_INPUT_E << " approximate energies (in keV) corresponding to peaks used to the fit." << endl;
        cout << "Default values will be used if arguments (other than analysis_tree) are omitted." << endl;
        return 0;
    }else if((argc >= 8)&&(argc < (5+MAX_INPUT_E))){
        sfile = argv[1];
        evalWindow = atof(argv[2]);
        tWindow = atof(argv[3]);
        soutfile = argv[4];
        for(int arg=5; arg<argc; arg++){
            corrEnergy[arg-5] = atof(argv[arg]);
        }
        numEnVals = argc-5;
    }else{
        printf("Incorrect arguments\nArguments: ReprocessGain_SMOL smol_file initial_time_window corr_time_window energy1 energy2 energy3 output_smolfile_suffix\n");
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

    cout << "Gain evaluation time window: " << evalWindow << " sec" << endl;
    cout << "Gain correction time window: " << tWindow << " sec" << endl;
    printf("Gain correction energies: ");
    for(int i=0;i<(numEnVals-1);i++){
        printf("%0.3f, ",corrEnergy[i]);
    }
    printf("%0.3f keV\n",corrEnergy[numEnVals-1]);
    cout << "Output file: " << soutfile << endl;

    memset(actualEnergy,0,sizeof(actualEnergy)); //zero out
    actualEnDetermined = 0;

    const char *dot = strrchr(sfile, '.'); //get the file extension
    if(dot==NULL){
        cout << "ERROR: couldn't get SMOL tree or list file name." << endl;
        return 0;
    }

    if(strcmp(dot + 1, "smole6") == 0){
        snprintf(outName,63,"%s_%s.smole6",sfile,soutfile);
        mysort->SortData(sfile, outName, evalWindow, tWindow, corrEnergy, numEnVals);
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
                    snprintf(outName,63,"%s_%s.smole6",str,soutfile);
                    mysort->SortData(str, outName, evalWindow, tWindow, corrEnergy, numEnVals);
                }
            }
        }
    }else{
        cout << "ERROR: improper file extension for SMOL tree or list (should be .smole6 or .list)." << endl;
        return 0;
    }

    

    return 0;
}
