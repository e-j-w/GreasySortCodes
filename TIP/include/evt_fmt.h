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

#define MAX_EVT_HIT 16

typedef struct
{
  uint8_t numTigHits;
	uint8_t numNoABHits;
	uint8_t numCsIHits;
	uint8_t numRFHits;
	double evtTimeNs; //time of the first hit, in ns - float doesn't have enough precision
}evt_header;

typedef struct
{
	float timeOffsetNs; //relative to evtTimeNs
	float energy;
	uint8_t core;
	uint8_t seg;
}tig_hit;

typedef struct
{
	float timeOffsetNs; //relative to evtTimeNs
	float energy;
	float PID;
	uint8_t detNum;
	uint8_t metadata; //bits 0-2: CsI fit type
}csi_hit;

typedef struct
{
	float timeOffsetNs; //relative to evtTimeNs
	float freq;
	float phase;
}rf_hit;

//struct to hold an event resident in memory,
//not used in the actual file written to disk
typedef struct
{
  evt_header header;
	tig_hit tigHit[MAX_EVT_HIT]; //addback
	tig_hit noABHit[MAX_EVT_HIT]; //non-addback
	csi_hit csiHit[MAX_EVT_HIT];
	rf_hit rfHit;
}sorted_evt;

inline static double tigHitTime(const sorted_evt *evt, const int hit){return (double)(evt->header.evtTimeNs + evt->tigHit[hit].timeOffsetNs);};
inline static double noABHitTime(const sorted_evt *evt, const int hit){return (double)(evt->header.evtTimeNs + evt->noABHit[hit].timeOffsetNs);};
inline static double csiHitTime(const sorted_evt *evt, const int hit){return (double)(evt->header.evtTimeNs + evt->csiHit[hit].timeOffsetNs);};
inline static uint8_t csiFitType(const csi_hit *hit){return (uint8_t)(hit->metadata & 7U);};

//reads an event into the sortedEvt struct,
//returns the number of events read
static int readSMOLEvent(FILE *inp, sorted_evt *sortedEvt){
	//read event
	memset(sortedEvt,0,sizeof(sorted_evt));
	fread(&sortedEvt->header,sizeof(evt_header),1,inp);
	for(int i = 0; i<sortedEvt->header.numTigHits;i++){
		fread(&sortedEvt->tigHit[i],sizeof(tig_hit),1,inp);
	}
	for(int i = 0; i<sortedEvt->header.numNoABHits;i++){
		fread(&sortedEvt->noABHit[i],sizeof(tig_hit),1,inp);
	}
	for(int i = 0; i<sortedEvt->header.numCsIHits;i++){
		fread(&sortedEvt->csiHit[i],sizeof(csi_hit),1,inp);
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
