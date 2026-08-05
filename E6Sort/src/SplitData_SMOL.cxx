//Sort code to remove crystals from SMOL files

#define SplitData_SMOL_cxx
#include "common.h"
#include "lin_eq_solver.h"
#include "evt_fmt.h"
#include "SplitData_SMOL.h"

using namespace std;

FILE *inp, *out;
char outName[512];

void SortData(const char *sfile, const char *outfile, const int numSplit){

    uint64_t sentries = 0U;
    uint64_t pileupCtrs[16];
    sorted_evt sortedEvt;
    uint64_t hitBuildFlags = 0;
    uint64_t numSeparatedEvents = 0;

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
    memset(pileupCtrs,0,sizeof(pileupCtrs)); //not used in the output files, zero out to avoid confusion later
    sentries &= 0xFFFFFFFFFFFF; // only first 48 bits specify number of events

    out = NULL;
    int outFileInd = 0;
    const uint64_t entriesPerFile = (uint64_t)(ceil((sentries*1.0)/(1.0*numSplit)));

    uint64_t outSepEntries = 0;
    for(Long64_t jentry = 0; jentry < sentries; jentry++){

        //setup the output file if needed
        if(out == NULL){
            snprintf(outName,511,"%s_%i.smol",outfile,outFileInd);
            out = fopen(outName, "wb");
            if(out == NULL){
                cout << "ERROR: couldn't open output file " << outName << endl;
                return;
            }
            printf("Will write %lu events to file: %s\n",entriesPerFile,outName);
            fwrite(&entriesPerFile,sizeof(uint64_t),1,out);
            if(smolVersion > 0){
                fwrite(&pileupCtrs,sizeof(pileupCtrs),1,out);
            }
        }

        //read event from input file
        if(readSMOLEvent(inp,&sortedEvt)==0){
            cout << "ERROR: bad event data in entry " << jentry << "." << endl;
            exit(-1);
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
        outSepEntries++;

        if(outSepEntries >= entriesPerFile){
            //write the number of separated events to the beginning of the file
            fseek(out,0,SEEK_SET);
            outSepEntries &= 0xFFFFFFFFFFFF; // only first 48 bits specify number of events
            uint64_t smolFormatVersion = 1;
            outSepEntries |= (smolFormatVersion << 48);
            fwrite(&outSepEntries,sizeof(uint64_t),1,out);
            printf(" Wrote %lu separated events to: %s\n",outSepEntries & 0xFFFFFFFFFFFF,outName);
            fclose(out);
            out = NULL; //to make sure a new file is opened during the next loop iteration
            outFileInd++;
            outSepEntries = 0;
        }

        if(outFileInd > numSplit){
            printf("ERROR: expected number of output files exceeded!\n"); //regulations are written in blood...
            return;
        }
        

        if(jentry % 9713 == 0)
            cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;

    }
    cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;

    if(out != NULL){
        //write the number of separated events to the beginning of the file
        fseek(out,0,SEEK_SET);
        outSepEntries &= 0xFFFFFFFFFFFF; // only first 48 bits specify number of events
        uint64_t smolFormatVersion = 1;
        outSepEntries |= (smolFormatVersion << 48);
        fwrite(&outSepEntries,sizeof(uint64_t),1,out);
        printf("Wrote %lu separated events to: %s\n",outSepEntries & 0xFFFFFFFFFFFF,outName);
        fclose(out);
    }
    

    fclose(inp);

    return;

}
int main(int argc, char **argv){

    char const *sfile;
    char const *soutfile;
    char outName[256];
    int numSplit = 0;

    if(argc <= 1){
        cout << "Arguments: SplitData_SMOL smol_file output_smolfile_suffix num_split" << endl;
        cout << "A code for splitting SMOL trees into sub-trees (eg. to analyze data at a sub-run level)." << endl;
        cout << "  *smol_file* is a single SMOL tree (extension .smol)." << endl;
        cout << "  *num_split* is the number of (equally sized) trees to split the input file into." << endl;
        return 0;
    }else if(argc == 4){
        sfile = argv[1];
        soutfile = argv[2];
        numSplit = atoi(argv[3]);
    }else{
        printf("Incorrect arguments\nArguments: SplitData_SMOL smol_file output_smolfile_suffix num_split\n");
        return 0;
    }

    printf("Starting SplitData_SMOL code\n");

    if(strcmp(soutfile,"")==0){
        cout << "ERROR: output suffix cannot be empty." << endl;
        return 0;
    }
    cout << "Output file suffix: " << soutfile << endl;

    char filePrefix[256];

    const char *dot = strrchr(sfile, '.'); //get the file extension
    if(dot==NULL){
        cout << "ERROR: couldn't get SMOL tree file name." << endl;
        return 0;
    }

    if(numSplit < 2){
        printf("ERROR: must split the input file into at least 2 files.\n");
        return 0;
    }
    printf("Will split tree into %i files.\n",numSplit);

    if(strcmp(dot + 1, "smol") == 0){
        strncpy(filePrefix,sfile,256);
        const char *tok = strtok(filePrefix,"."); //get the filename without the extension
        if(tok!=NULL){
            snprintf(outName,255,"%s_%s",basename(tok),soutfile);
            //printf("Will write to file: %s\n",outName);
            SortData(sfile, outName, numSplit);
        }else{
            cout << "ERROR: improperly formatted filename: " << sfile << endl;
            return 0;
        }
    }else{
        cout << "ERROR: improper file extension for *smol_file* argument (should be .smol)." << endl;
        return 0;
    }

    return 0;
}
