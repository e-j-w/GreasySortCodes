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
  uint8_t numTigABHits;
	uint8_t numCsIHits;
	uint8_t numRFHits;
}evt_header;

typedef struct
{
	double timeNs; //absolute time - float doesn't have enough precision
	float energy;
	uint8_t core;
	uint8_t seg;
}tigab_hit;

typedef struct
{
	double timeNs; //absolute time - float doesn't have enough precision
	float energy;
	float PID;
	uint8_t detNum;
}csi_hit;

typedef struct
{
	double timeNs; //absolute time - float doesn't have enough precision
	float freq;
	float phase;
}rf_hit;

//struct to hold an event resident in memory,
//not used in the actual file written to disk
typedef struct
{
  evt_header header;
	tigab_hit tigHit[MAX_EVT_HIT];
	csi_hit csiHit[MAX_EVT_HIT];
	rf_hit rfHit;
}sorted_evt;

#endif
