#ifndef evtfmt_h
#define evtfmt_h

//structs defining the intermediate tree format used for sorting GRIFFIN data
//
//see smol-format.h in midas2smol

#include <stdint.h> //allows uint8_t and similiar types

#define MAX_EVT_HIT 64

typedef struct
{
	uint8_t numNoABHits;
	uint8_t metadata; //bit 0: TIGRESS (on) or GRIFFIN (off), bit 1: any suppressor fired, bit 7: always set, for data validation
	double evtTimeNs; //time of the first hit, in ns - float doesn't have enough precision
}evt_header;

typedef struct
{
	float timeOffsetNs; //relative to evtTimeNs
	float energy;
	uint8_t tsDiff;
	uint8_t core; //0-indexed
}hpge_hit;

//struct to hold an event resident in memory,
//not used in the actual file written to disk
typedef struct
{
	evt_header header;
	hpge_hit noABHit[MAX_EVT_HIT]; //non-addback
}sorted_evt;

inline static double noABHitTime(const sorted_evt *evt, const int hit){return (double)(evt->header.evtTimeNs + evt->noABHit[hit].timeOffsetNs);};
inline static double isTIGRESSData(const sorted_evt *evt){return evt->header.metadata & 1U;};

//reads an event into the sortedEvt struct,
//returns the number of events read
static int readSMOLEvent(FILE *inp, sorted_evt *sortedEvt){
	//read event
	memset(sortedEvt,0,sizeof(sorted_evt));
	fread(&sortedEvt->header,sizeof(evt_header),1,inp);
	for(uint8_t i = 0; i<sortedEvt->header.numNoABHits;i++){
		fread(&sortedEvt->noABHit[i].timeOffsetNs,sizeof(float),1,inp);
		fread(&sortedEvt->noABHit[i].energy,sizeof(float),1,inp);
		fread(&sortedEvt->noABHit[i].tsDiff,sizeof(uint8_t),1,inp);
		fread(&sortedEvt->noABHit[i].core,sizeof(uint8_t),1,inp);
	}
	if(!(sortedEvt->header.metadata & (1U << 7))){
		//validation bit not set, might be older format
		//that used a footer byte for validation
		uint8_t footerVal = 0;
		fread(&footerVal,sizeof(uint8_t),1,inp);
		if(footerVal != 227U){
			//invalid footer
			//dump bad event data
			printf("Bad event:\n");
			printf("  metadata: %u\n",sortedEvt->header.metadata);
			printf("  event time: %f\n",sortedEvt->header.evtTimeNs);
			printf("  num hits: %u\n",sortedEvt->header.numNoABHits);
			for(uint8_t i = 0; i<sortedEvt->header.numNoABHits;i++){
				printf("    hit %u - time offset: %f, energy: %f, core: %u\n",i,sortedEvt->noABHit[i].timeOffsetNs,sortedEvt->noABHit[i].energy,sortedEvt->noABHit[i].core);
			}
			return 0;
		}
	}
	return 1;
}

#endif
