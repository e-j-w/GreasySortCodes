#ifndef evtfmt_h
#define evtfmt_h

//structs defining the intermediate tree format used for sorting TIP data
//
//Data written to disk has the format:
//-64-bit number containing number of events
//-Individual event data, containing:
//  -evt_header struct
//  -one tigab_hit struct for each addback hit
//  -one csi_hit struct for each csi hit
//  -one rf_hit struct for each rf hit
//  -footer containing magic uint8_t value (227U), used for data validation

#include <stdint.h> //allows uint8_t and similiar types

#define MAX_EVT_HIT 64

typedef struct
{
	uint8_t numABHits;
	uint8_t numNoABHits;
	uint8_t metadata; //bit 0: TIGRESS (on) or GRIFFIN (off), bit 1: any suppressor fired
	double evtTimeNs; //time of the first hit, in ns - float doesn't have enough precision
}evt_header;

typedef struct
{
	float timeOffsetNs; //relative to evtTimeNs
	float energy;
	uint8_t core; //0-indexed
}hpge_hit;

//struct to hold an event resident in memory,
//not used in the actual file written to disk
typedef struct
{
	evt_header header;
	hpge_hit ABHit[MAX_EVT_HIT]; //addback
	hpge_hit noABHit[MAX_EVT_HIT]; //non-addback
}sorted_evt;

inline static double ABHitTime(const sorted_evt *evt, const int hit){return (double)(evt->header.evtTimeNs + evt->ABHit[hit].timeOffsetNs);};
inline static double noABHitTime(const sorted_evt *evt, const int hit){return (double)(evt->header.evtTimeNs + evt->noABHit[hit].timeOffsetNs);};
inline static double isTIGRESSData(const sorted_evt *evt){return evt->header.metadata & 1U;};

//reads an event into the sortedEvt struct,
//returns the number of events read
static int readSMOLEvent(FILE *inp, sorted_evt *sortedEvt){
	//read event
	memset(sortedEvt,0,sizeof(sorted_evt));
	fread(&sortedEvt->header,sizeof(evt_header),1,inp);
	for(int i = 0; i<sortedEvt->header.numABHits;i++){
		fread(&sortedEvt->ABHit[i].timeOffsetNs,sizeof(float),1,inp);
		fread(&sortedEvt->ABHit[i].energy,sizeof(float),1,inp);
		fread(&sortedEvt->ABHit[i].core,sizeof(uint8_t),1,inp);
	}
	for(int i = 0; i<sortedEvt->header.numNoABHits;i++){
		fread(&sortedEvt->noABHit[i].timeOffsetNs,sizeof(float),1,inp);
		fread(&sortedEvt->noABHit[i].energy,sizeof(float),1,inp);
		fread(&sortedEvt->noABHit[i].core,sizeof(uint8_t),1,inp);
	}
	uint8_t footerVal = 0;
	fread(&footerVal,sizeof(uint8_t),1,inp);
	if(footerVal != 227U){
		//invalid footer
		return 0;
	}
	return 1;
}

#endif
