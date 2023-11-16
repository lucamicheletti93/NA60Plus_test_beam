#!/bin/bash 
cp CustomModules/ClusteringSpatial.cpp /opt/corryvreckan/src/modules/ClusteringSpatial/ClusteringSpatial.cpp
cp CustomModules/ClusteringSpatial.h /opt/corryvreckan/src/modules/ClusteringSpatial/ClusteringSpatial.h

cp CustomModules/Tracking4D.cpp /opt/corryvreckan/src/modules/Tracking4D/Tracking4D.cpp
cp CustomModules/Tracking4D.h /opt/corryvreckan/src/modules/Tracking4D/Tracking4D.h

cp CustomModules/Track.hpp /opt/corryvreckan/src/objects/Track.hpp

cp CustomModules/Linkdef.h /opt/corryvreckan/src/objects/Linkdef.h
cp CustomModules/CorryTree.hpp /opt/corryvreckan/src/objects/CorryTree.hpp
cp CustomModules/CMakeLists.txt /opt/corryvreckan/src/objects/CMakeLists.txt

cd /opt/corryvreckan/build/ && make install -j12 && cd /local
