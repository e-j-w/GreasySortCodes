//File contains common functions which all other sort codes
//have access to

#define common_cxx
#include "../include/common.h"
#include "position_arrays.cxx"

using namespace std;

TVector3 getTigVector(uint8_t core, uint8_t seg){
  TVector3 hitPos(0,0,0);
  if(core > 63){
    printf("WARNING: bad core value (%u)\n",core);
  }else if(seg > 8){
    printf("WARNING: bad segment value (%u)\n",seg);
  }else{
    switch(core % 4){
      case 3:
        //white core
        hitPos.SetXYZ(GeWhitePositionBack[(core/4) + 1][seg][0],GeWhitePositionBack[(core/4) + 1][seg][1],GeWhitePositionBack[(core/4) + 1][seg][2]);
        break;
      case 2:
        //red core
        hitPos.SetXYZ(GeRedPositionBack[(core/4) + 1][seg][0],GeRedPositionBack[(core/4) + 1][seg][1],GeRedPositionBack[(core/4) + 1][seg][2]);
        break;
      case 1:
        //green core
        hitPos.SetXYZ(GeGreenPositionBack[(core/4) + 1][seg][0],GeGreenPositionBack[(core/4) + 1][seg][1],GeGreenPositionBack[(core/4) + 1][seg][2]);
        break;
      case 0:
      default:
        //blue core
        hitPos.SetXYZ(GeBluePositionBack[(core/4) + 1][seg][0],GeBluePositionBack[(core/4) + 1][seg][1],GeBluePositionBack[(core/4) + 1][seg][2]);
        break;
    }
  }
  return hitPos;
}

bool gate1D(const double value, const double min, const double max){
  if (min < value && value < max)
    return true;
  else
    return false;
}


int getTIGRESSRing(const float theta){
  if((theta < 0.0f)||(theta > 180.0f)){
    cout << "Invalid TIGRESS angle: " << theta << endl;
    return 0;
  }else if(theta <= 45.0f){
    return 0;
  }else if(theta <= 66.0f){
    return 1;
  }else if(theta <= 90.0f){
    return 2;
  }else if(theta <= 112.0f){
    return 3;
  }else if(theta <= 135.0f){
    return 4;
  }else{
    return 5;
  }
}


//Get the segment 'ring'
//there are actually 8 unique segment angles per core, clustered
//over 'rings' each with 4 angles covering a ~2.5 deg range
//there are also core-only hits without any assigned segment
int getTIGRESSSegmentRing(const float theta){
  if((theta >= 30.0f)&&(theta <= 37.0f)){
    return 0; //34.325 deg (low 32.95, high 35.55)
  }else if((theta >= 40.0f)&&(theta <= 43.0f)){
    return 1; //41.283 deg (low 40.55, high 42.15, missing 1 angle)
  }else if((theta >= 48.0f)&&(theta <= 52.0f)){
    return 2; //49.875 deg (low 49.15, high 50.55)
  }else if((theta >= 55.0f)&&(theta <= 59.0f)){
    return 3; //57.075 deg (low 56.05, high 58.25)
  }else if((theta >= 77.0f)&&(theta <= 80.0f)){
    return 4; //78.50 deg (low 77.75, high 79.35)
  }else if((theta >= 84.0f)&&(theta <= 88.0f)){
    return 5; //85.683 deg (low 85.35, high 86.05, missing 1 angle)
  }else if((theta >= 92.0f)&&(theta <= 96.0f)){
    return 6; //94.317 deg (low 93.95, high 94.65, missing 1 angle)
  }else if((theta >= 100.0f)&&(theta <= 103.0f)){
    return 7; //101.50 deg (low 100.65, high 102.25)
  }else if((theta >= 121.0f)&&(theta <= 125.0f)){
    return 8; //122.920 deg (low 121.75, high 123.95)
  }else if((theta >= 128.0f)&&(theta <= 132.0f)){
    return 9; //130.120 deg (low 129.45, high 130.85)
  }else if((theta >= 137.0f)&&(theta <= 141.0f)){
    return 10; //138.72 deg (low 137.85, high 139.45, missing 1 angle)
  }else if((theta >= 143.0f)&&(theta <= 148.0f)){
    return 11; //145.680 deg (low 144.45, high 147.05)
  }else{
    //cout << "Invalid TIGRESS angle: " << theta << endl;
    return 12; //no segment
  }
}

int getTIGRESSPhiRing(const float theta, const float phi){
  int ring = getTIGRESSSegmentRing(theta);
  if(ring < 12){
    if((phi > -180)&&(phi < -167)){
      return 0; 
    }else if((phi > -165)&&(phi < -160)){
      return 1; 
    }else if((phi > -155)&&(phi < -148)){
      return 2; 
    }else if((phi > -148)&&(phi < -134)){
      return 3;
    }else if((phi > -127)&&(phi < -122)){
      return 4;
    }else if((phi > -118)&&(phi < -114)){
      return 5;
    }else if((phi > -110)&&(phi < -106)){
      return 6;
    }else if((phi > -103)&&(phi < -100)){
      return 7;
    }else if((phi > -91)&&(phi < -77)){
      return 8;
    }else if((phi > -77)&&(phi < -70)){
      return 9;
    }else if((phi > -65)&&(phi < -58)){
      return 10;
    }else if((phi > -58)&&(phi < -45)){
      return 11;
    }else if((phi > -36)&&(phi < -32)){
      return 12;
    }else if((phi > -28)&&(phi < -25)){
      return 13;
    }else if((phi > -20)&&(phi < -16)){
      return 14;
    }else if((phi > -13)&&(phi < -8)){
      return 15;
    }else if((phi > 0)&&(phi < 13)){
      return 16;
    }else if((phi > 13)&&(phi < 22)){
      return 17;
    }else if((phi > 22)&&(phi < 32)){
      return 18;
    }else if((phi > 32)&&(phi < 46)){
      return 19;
    }else if((phi > 53)&&(phi < 59)){
      return 20;
    }else if((phi > 61)&&(phi < 66)){
      return 21;
    }else if((phi > 69)&&(phi < 74)){
      return 22;
    }else if((phi > 77)&&(phi < 82)){
      return 23;
    }else if((phi > 89)&&(phi < 103)){
      return 24;
    }else if((phi > 103)&&(phi < 114)){
      return 25;
    }else if((phi > 114)&&(phi < 122)){
      return 26;
    }else if((phi > 122)&&(phi < 136)){
      return 27;
    }else if((phi > 143)&&(phi < 149)){ //estimated, no detector present here during S2182
      return 28;
    }else if((phi > 151)&&(phi < 156)){ //estimated, no detector present here during S2182
      return 29;
    }else if((phi > 159)&&(phi < 164)){ //estimated, no detector present here during S2182
      return 30;
    }else if((phi > 167)&&(phi < 172)){ //estimated, no detector present here during S2182
      return 31;
    }else{
      //cout << "Invalid TIGRESS angle: " << theta << endl;
      return 32; //no phi ring
    }
  }else{
    //cout << "Invalid TIGRESS angle: " << theta << endl;
    return 32; //no phi ring
  }
}

