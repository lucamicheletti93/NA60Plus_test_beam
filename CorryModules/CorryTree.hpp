#include <vector>

class CorryTrack {
    public:
        float vx; // vector x component 
        float vy; // vector y component 
        //float vz; // vector z component not needed -> always 1 
        float px; //x position
        float py; //y position
        //float pz; //z position -> always 1
        float chi2ndof; 
        //Clusters info
        //ALPIDE 1
        int   size1; //cluster size
        //ALPIDE 2
        int   size2; //cluster size
        //ALPIDE 3
        int   size3; //cluster size
        //ALPIDE 4
        int   size4; //cluster size
        //ALPIDE 5
        int   size5; //cluster size
        //ALPIDE 6
        int   size6; //cluster size
};

class CorryBeam {
    public:
        float vx; // vector x component 
        float vy; // vector y component 
        float px; //x position
        float py; //y position
        int size1; //cluster size of the largest cluster
        int size2; //cluster size of the largest cluster
};

class CorryEvent : public TObject {
    public:
        CorryBeam beam;
        std::vector<CorryTrack> tracks;
    ClassDef(CorryEvent, 1)
};

ClassImp(CorryEvent)