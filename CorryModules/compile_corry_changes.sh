#!/bin/bash 
cp ClusteringSpatial.cpp /opt/corryvreckan/src/modules/ClusteringSpatial/ClusteringSpatial.cpp
cp ClusteringSpatial.h /opt/corryvreckan/src/modules/ClusteringSpatial/ClusteringSpatial.h

cp Tracking4D.cpp /opt/corryvreckan/src/modules/Tracking4D/Tracking4D.cpp
cp Tracking4D.h /opt/corryvreckan/src/modules/Tracking4D/Tracking4D.h

cd /opt/corryvreckan/build/ && make install -j12 && cd /local