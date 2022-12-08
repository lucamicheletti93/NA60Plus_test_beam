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
        float clx1; //x position
        float cly1; //y position
        //ALPIDE 2
        int   size2; //cluster size
        float clx2; //x position
        float cly2; //y position
        //ALPIDE 3
        int   size3; //cluster size
        float clx3; //x position
        float cly3; //y position
        //ALPIDE 4
        int   size4; //cluster size
        float clx4; //x position
        float cly4; //y position
        //ALPIDE 5
        int   size5; //cluster size
        float clx5; //x position
        float cly5; //y position
};

class CorryBeam {
    public:
        float x; //x position of the largest cluster
        float y; //y position of the largest cluster
        int size; //cluster size of the largest cluster
};

class CorryEvent : public TObject {
    public:
        CorryBeam beam;
        std::vector<CorryTrack> tracks;
    ClassDef(CorryEvent, 1)
};

ClassImp(CorryEvent)