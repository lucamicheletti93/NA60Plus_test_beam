/**
 * @file
 * @brief Implementation of module Tracking4D
 *
 * @copyright Copyright (c) 2017-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "Tracking4D.h"
#include <TCanvas.h>
#include <TDirectory.h>

#include "tools/cuts.h"
#include "tools/kdtree.h"

using namespace corryvreckan;
using namespace std;

Tracking4D::Tracking4D(Configuration& config, std::vector<std::shared_ptr<Detector>> detectors)
    : Module(config, std::move(detectors)) {

    // Backwards compatibility: also allow timing_cut to be used for time_cut_abs and spatial_cut for spatial_cut_abs
    config_.setAlias("spatial_cut_abs", "spatial_cut", true);

    config_.setDefault<size_t>("min_hits_on_track", 6);
    config_.setDefault<bool>("exclude_dut", true);
    config_.setDefault<bool>("apply_selections", false);
    config_.setDefault<double>("max_plot_chi2", 50.0);
    config_.setDefault<bool>("reject_by_roi", false);
    config_.setDefault<bool>("unique_cluster_usage", false);
    config_.setDefault<float>("max_delta_eta", 1000);
    config_.setDefault<float>("min_eta_y", 0);
    config_.setDefault<float>("min_eta_x", 0);
    config_.setDefault<float>("max_eta", 1000);
    config_.setDefault<float>("target_position", Units::get<float>(-1, "cm"));
    config_.setDefault<int>("max_range_multiplicity", 100);
    config_.setDefault<float>("spatial_cuts_max", 1000);
    
    if(config_.count({"spatial_cut_rel", "spatial_cut_abs"}) == 0) {
        config_.setDefault("spatial_cut_rel", 3.0);
    }


    spatial_cuts_tmp_ = corryvreckan::calculate_cut<XYVector>("spatial_cut", config_, get_regular_detectors(true));


    min_hits_on_track_ = config_.get<size_t>("min_hits_on_track");
    exclude_DUT_ = config_.get<bool>("exclude_dut");

    require_detectors_ = config_.getArray<std::string>("require_detectors", {});
    exclude_from_seed_ = config_.getArray<std::string>("exclude_from_seed", {});
    max_plot_chi2_ = config_.get<double>("max_plot_chi2");
    reject_by_ROI_ = config_.get<bool>("reject_by_roi");
    unique_cluster_usage_ = config_.get<bool>("unique_cluster_usage");
    apply_selections_ = config_.get<bool>("apply_selections");
    max_delta_eta_ = config_.get<float>("max_delta_eta");
    min_eta_y_ = config_.get<float>("min_eta_y");
    min_eta_x_ = config_.get<float>("min_eta_x");
    max_eta_ = config_.get<float>("max_eta");
    spatial_cuts_max_ = config_.get<float>("spatial_cuts_max");
    target_position_ = config_.get<float>("target_position");
    max_range_multiplicity_ = config_.get<int>("max_range_multiplicity");
}

void Tracking4D::initialize() {

    // Set up histograms
    std::string title = "Track #chi^{2};#chi^{2};events";
    trackChi2 = new TH1F("trackChi2", title.c_str(), 300, 0, 3 * max_plot_chi2_);
    title = "Track #chi^{2}/ndof;#chi^{2}/ndof;events";
    trackChi2ndof = new TH1F("trackChi2ndof", title.c_str(), 500, 0, max_plot_chi2_);
    title = "Clusters per track;clusters;tracks";
    clustersPerTrack = new TH1F("clustersPerTrack", title.c_str(), 10, -0.5, 9.5);
    title = "Track multiplicity;tracks;events";
    tracksPerEvent = new TH1F("tracksPerEvent", title.c_str(), max_range_multiplicity_, -0.5, max_range_multiplicity_-0.5);
    title = "Track angle X;angle_{x} [rad];events";
    trackAngleX = new TH1F("trackAngleX", title.c_str(), 2000, -3.14/2, 3.14/2);
    title = "Track angle Y;angle_{y} [rad];events";
    trackAngleY = new TH1F("trackAngleY", title.c_str(), 2000, -3.14/2, 3.14/2);


    hFilter = new TH1F("FilteredEvents", "Events filtered;events", 3, 0.5, 3.5);
    hFilter->GetXaxis()->SetBinLabel(1, "delta eta");
    hFilter->GetXaxis()->SetBinLabel(2, "eta min");
    hFilter->GetXaxis()->SetBinLabel(3, "eta max");

    event_info = new TTree("tree","event info");
    event_corry = new CorryEvent;
    event_info->Branch("event",&event_corry);
    int plane_counter = 0;
    // Loop over all planes
    for(auto& detector : get_regular_detectors(true)) {
        auto detectorID = detector->getName();
        
        if(plane_counter < 3)
            spatial_cuts_[detectorID] = static_cast<float>(spatial_cuts_tmp_[detector].x());
        else
            spatial_cuts_[detectorID] = static_cast<float>(spatial_cuts_tmp_[detector].x()*(plane_counter-3.)/4.);//hardcoded stuff
        plane_counter++;
        
        TDirectory* directory = getROOTDirectory();
        TDirectory* local_directory = directory->mkdir(detectorID.c_str());

        if(local_directory == nullptr) {
            throw RuntimeError("Cannot create or access local ROOT directory for module " + this->getUniqueName());
        }
        local_directory->cd();

        local_intersects_[detectorID] = new TH2F("local_intersect",
                                                 "local intersect, col, row",
                                                 detector->nPixels().X(),
                                                 0,
                                                 detector->nPixels().X(),
                                                 detector->nPixels().Y(),
                                                 0,
                                                 detector->nPixels().Y());
        global_intersects_[detectorID] = new TH2F("global_intersect",
                                                  "global intersect, global intercept x [mm];global intercept y [mm]",
                                                  600,
                                                  -30,
                                                  30,
                                                  600,
                                                  -30,
                                                  30);

        // Do not create plots for detectors not participating in the tracking:
        if(exclude_DUT_ && detector->isDUT()) {
            continue;
        }
        // local
        TDirectory* local_res = local_directory->mkdir("local_residuals");
        local_res->cd();
        title = detectorID + "Local Residual X;x-x_{track} [mm];events";
        residualsX_local[detectorID] =
            new TH1F("LocalResidualsX", title.c_str(), 500, -3 * detector->getPitch().X(), 3 * detector->getPitch().X());
        title = detectorID + "Local  Residual X, cluster column width 1;x-x_{track} [mm];events";
        residualsXwidth1_local[detectorID] = new TH1F(
            "LocalResidualsXwidth1", title.c_str(), 500, -3 * detector->getPitch().X(), 3 * detector->getPitch().X());
        title = detectorID + "Local  Residual X, cluster column width  2;x-x_{track} [mm];events";
        residualsXwidth2_local[detectorID] = new TH1F(
            "LocalResidualsXwidth2", title.c_str(), 500, -3 * detector->getPitch().X(), 3 * detector->getPitch().X());
        title = detectorID + "Local  Residual X, cluster column width  3;x-x_{track} [mm];events";
        residualsXwidth3_local[detectorID] = new TH1F(
            "LocalResidualsXwidth3", title.c_str(), 500, -3 * detector->getPitch().X(), 3 * detector->getPitch().X());
        title = detectorID + "Local  Residual Y;y-y_{track} [mm];events";
        residualsY_local[detectorID] =
            new TH1F("LocalResidualsY", title.c_str(), 500, -3 * detector->getPitch().Y(), 3 * detector->getPitch().Y());
        title = detectorID + "Local  Residual Y, cluster row width 1;y-y_{track} [mm];events";
        residualsYwidth1_local[detectorID] = new TH1F(
            "LocalResidualsYwidth1", title.c_str(), 500, -3 * detector->getPitch().Y(), 3 * detector->getPitch().Y());
        title = detectorID + "Local  Residual Y, cluster row width 2;y-y_{track} [mm];events";
        residualsYwidth2_local[detectorID] = new TH1F(
            "LocalResidualsYwidth2", title.c_str(), 500, -3 * detector->getPitch().Y(), 3 * detector->getPitch().Y());
        title = detectorID + "Local  Residual Y, cluster row width 3;y-y_{track} [mm];events";
        residualsYwidth3_local[detectorID] = new TH1F(
            "LocalResidualsYwidth3", title.c_str(), 500, -3 * detector->getPitch().Y(), 3 * detector->getPitch().Y());

        title = detectorID + " Pull X;x-x_{track}/resolution;events";
        pullX_local[detectorID] = new TH1F("LocalpullX", title.c_str(), 500, -5, 5);

        title = detectorID + " Pull Y;y-y_{track}/resolution;events";
        pullY_local[detectorID] = new TH1F("Localpully", title.c_str(), 500, -5, 5);
        // global
        TDirectory* global_res = local_directory->mkdir("global_residuals");
        global_res->cd();
        title = detectorID + "global Residual X;x-x_{track} [mm];events";
        residualsX_global[detectorID] =
            new TH1F("GlobalResidualsX", title.c_str(), 500, -3 * detector->getPitch().X(), 3 * detector->getPitch().X());

        title = detectorID + " global  Residual X vs. global position X;x-x_{track} [mm];x [mm]";
        residualsX_vs_positionX_global[detectorID] = new TH2F("GlobalResidualsX_vs_GlobalPositionX",
                                                              title.c_str(),
                                                              500,
                                                              -3 * detector->getPitch().X(),
                                                              3 * detector->getPitch().X(),
                                                              400,
                                                              -detector->getSize().X() / 1.5,
                                                              detector->getSize().X() / 1.5);
        title = detectorID + " global  Residual X vs. global position Y;x-x_{track} [mm];y [mm]";
        residualsX_vs_positionY_global[detectorID] = new TH2F("GlobalResidualsX_vs_GlobalPositionY",
                                                              title.c_str(),
                                                              500,
                                                              -3 * detector->getPitch().X(),
                                                              3 * detector->getPitch().X(),
                                                              400,
                                                              -detector->getSize().Y() / 1.5,
                                                              detector->getSize().Y() / 1.5);

        title = detectorID + "global  Residual X, cluster column width 1;x-x_{track} [mm];events";
        residualsXwidth1_global[detectorID] = new TH1F(
            "GlobalResidualsXwidth1", title.c_str(), 500, -3 * detector->getPitch().X(), 3 * detector->getPitch().X());
        title = detectorID + "global  Residual X, cluster column width  2;x-x_{track} [mm];events";
        residualsXwidth2_global[detectorID] = new TH1F(
            "GlobalResidualsXwidth2", title.c_str(), 500, -3 * detector->getPitch().X(), 3 * detector->getPitch().X());
        title = detectorID + "global  Residual X, cluster column width  3;x-x_{track} [mm];events";
        residualsXwidth3_global[detectorID] = new TH1F(
            "GlobalResidualsXwidth3", title.c_str(), 500, -3 * detector->getPitch().X(), 3 * detector->getPitch().X());
        title = detectorID + " Pull X;x-x_{track}/resolution;events";
        pullX_global[detectorID] = new TH1F("GlobalpullX", title.c_str(), 500, -5, 5);
        title = detectorID + "global  Residual Y;y-y_{track} [mm];events";
        residualsY_global[detectorID] =
            new TH1F("GlobalResidualsY", title.c_str(), 500, -3 * detector->getPitch().Y(), 3 * detector->getPitch().Y());

        title = detectorID + " global  Residual Y vs. global position Y;y-y_{track} [mm];y [mm]";
        residualsY_vs_positionY_global[detectorID] = new TH2F("GlobalResidualsY_vs_GlobalPositionY",
                                                              title.c_str(),
                                                              500,
                                                              -3 * detector->getPitch().Y(),
                                                              3 * detector->getPitch().Y(),
                                                              400,
                                                              -detector->getSize().Y() / 1.5,
                                                              detector->getSize().Y() / 1.5);
        title = detectorID + " global  Residual Y vs. global position X;y-y_{track} [mm];x [mm]";
        residualsY_vs_positionX_global[detectorID] = new TH2F("GlobalResidualsY_vs_GlobalPositionX",
                                                              title.c_str(),
                                                              500,
                                                              -3 * detector->getPitch().Y(),
                                                              3 * detector->getPitch().Y(),
                                                              400,
                                                              -detector->getSize().X() / 1.5,
                                                              detector->getSize().X() / 1.5);

        title = detectorID + "global  Residual Y, cluster row width 1;y-y_{track} [mm];events";
        residualsYwidth1_global[detectorID] = new TH1F(
            "GlobalResidualsYwidth1", title.c_str(), 500, -3 * detector->getPitch().Y(), 3 * detector->getPitch().Y());
        title = detectorID + "global  Residual Y, cluster row width 2;y-y_{track} [mm];events";
        residualsYwidth2_global[detectorID] = new TH1F(
            "GlobalResidualsYwidth2", title.c_str(), 500, -3 * detector->getPitch().Y(), 3 * detector->getPitch().Y());
        title = detectorID + "global  Residual Y, cluster row width 3;y-y_{track} [mm];events";
        residualsYwidth3_global[detectorID] = new TH1F(
            "GlobalResidualsYwidth3", title.c_str(), 500, -3 * detector->getPitch().Y(), 3 * detector->getPitch().Y());
        title = detectorID + " Pull Y;y-y_{track}/resolution;events";
        pullY_global[detectorID] = new TH1F("Globalpully", title.c_str(), 500, -5, 5);;
    }
}
        
// Define your custom function
double Tracking4D::eta_threshold(double x, double params) {
    // Implement your function here
    double L = 2.5;
    double a = TMath::Tan(TMath::ASin((x * 2. + params) / (L * 2.)) / 2.);
    double b = TMath::Tan(TMath::ASin((x - params) / (L)) / 2.);
    return TMath::Abs(TMath::Log(a/b));
}

StatusCode Tracking4D::run(const std::shared_ptr<Clipboard>& clipboard) {

    //LOG(DEBUG) << "Start of event";
    // Container for all clusters, and detectors in tracking
    map<std::shared_ptr<Detector>, KDTree<Cluster>> trees;

    std::shared_ptr<Detector> reference_first, reference_last;
    for(auto& detector : get_regular_detectors(!exclude_DUT_)) {
        // Get the clusters
        auto tempClusters = clipboard->getData<Cluster>(detector->getName());
        //LOG(DEBUG) << "Detector " << detector->getName() << " has " << tempClusters.size() << " clusters on the clipboard";
        if(!tempClusters.empty()) {
            // Store them
            //LOG(DEBUG) << "Picked up " << tempClusters.size() << " clusters from " << detector->getName();

            trees.emplace(std::piecewise_construct, std::make_tuple(detector), std::make_tuple());
            trees[detector].buildTrees(tempClusters);

            // Get first and last detectors with clusters on them:
            if(std::find(exclude_from_seed_.begin(), exclude_from_seed_.end(), detector->getName()) ==
               exclude_from_seed_.end()) {
                if(!reference_first) {
                    reference_first = detector;
                }
                reference_last = detector;
            } else {
                //LOG(DEBUG) << "Not using " << detector->getName() << " as seed as chosen by config file.";
            }
        }
    }

    // If there are no detectors then stop trying to track
    if(trees.size() < 2) {
        // Fill histogram
        tracksPerEvent->Fill(0);

        //LOG(DEBUG) << "Too few hit detectors for finding a track; end of event.";
        return StatusCode::Success;
    }

    //LOG(DEBUG) << "Fitting the beam";
    // Output track container
    TrackVector tracks;
    CorryBeam beam;
    beam.vx = 0;
    beam.vy = 0;
    beam.px = 0;
    beam.py = 0;
    beam.size1 = 0;
    beam.size2 = 0;
    
    int beamMaxSize = 0;
    //LOG(DEBUG) << "Loop over ALPIDE 0";
    auto tempClusters = clipboard->getData<Cluster>("ALPIDE_0");
    //LOG(DEBUG) << "cluster abtained";
    Cluster* cluster1 = new Cluster();

    //LOG(DEBUG) << "cluster1 = temp[0]";
    for(auto tempCluster : tempClusters){
        beamMaxSize = static_cast<int>(tempCluster->size());
        if(beam.size1 < beamMaxSize){
            cluster1 = tempCluster.get();
            beam.size1 = beamMaxSize;
        }
    }

    //LOG(DEBUG) << "Loop over ALPIDE 1";
    tempClusters = clipboard->getData<Cluster>("ALPIDE_1");
    Cluster* cluster2 = new Cluster();
    for(auto tempCluster : tempClusters){
        beamMaxSize = static_cast<int>(tempCluster->size());
        if(beam.size2 < beamMaxSize){
            cluster2 = tempCluster.get();
            beam.size2 = beamMaxSize;
        }
    }

    //LOG(DEBUG) << "Add clusters "<< beam.size1 << " " << beam.size2;

    beam.vx = static_cast<float>((cluster2->global().x() - cluster1->global().x())/(cluster2->global().z() - cluster1->global().z()));
    beam.vy = static_cast<float>((cluster2->global().y() - cluster1->global().y())/(cluster2->global().z() - cluster1->global().z()));

    // Calculate the cluster (px, py) on the line
    beam.px = -beam.vx*static_cast<float>(cluster1->global().z())+static_cast<float>(cluster1->global().x());
    beam.py = -beam.vy*static_cast<float>(cluster1->global().z())+static_cast<float>(cluster1->global().y());

    event_corry->beam.vx = beam.vx;
    event_corry->beam.vy = beam.vy;
    event_corry->beam.px = beam.px;
    event_corry->beam.py = beam.py;
    event_corry->beam.size1 = beam.size1;
    event_corry->beam.size2 = beam.size2;

    float collision_x = beam.vx*100+beam.px;
    float collision_y = beam.vy*100+beam.px;

    //LOG(DEBUG) << "Beam fitted and saved";

    // Time cut for combinations of reference clusters and for reference track with additional detector
    for(auto& clusterFirst : trees[reference_first].getAllElements()) {

        for(auto& clusterLast : trees[reference_last].getAllElements()) {
            //LOG(DEBUG) << "Looking at next reference cluster pair";
            // The track finding is based on a straight line. Therefore a refTrack to extrapolate to the next plane is used
            StraightLineTrack refTrack;
            
            //measure eta

            double vx1 = (clusterFirst->global().x() - collision_x)/(clusterFirst->global().z() - target_position_);
            double vy1 = (clusterFirst->global().y() - collision_y)/(clusterFirst->global().z() - target_position_);

            double direction_tot = TMath::Sqrt(vx1*vx1+vy1*vy1+1*1);
            //double eta1 = TMath::ATanH(1/direction_tot); 

            double direction_x = TMath::Sqrt(vx1*vx1+1*1);
            double eta1_x = TMath::ATanH(1/direction_x); 

            double direction_y = TMath::Sqrt(vy1*vy1+1*1);
            double eta1_y = TMath::ATanH(1/direction_y); 

            double vx2 = (clusterLast->global().x() - collision_x)/(clusterLast->global().z() - target_position_);
            double vy2 = (clusterLast->global().y() - collision_y)/(clusterLast->global().z() - target_position_);

            direction_tot = TMath::Sqrt(vx2*vx2+vy2*vy2+1*1);
            //double eta2 = TMath::ATanH(1/direction_tot); //TMath::Log((direction_tot+1)/(direction_tot-1))/2.;


            direction_x = TMath::Sqrt(vx2*vx2+1*1);
            double eta2_x = TMath::ATanH(1/direction_x); 

            direction_y = TMath::Sqrt(vy2*vy2+1*1);
            double eta2_y = TMath::ATanH(1/direction_y); 
            //direction_x = TMath::Sqrt(vx2*vx2+1*1);
            //double eta2_x = TMath::ATanH(1/direction_x); 

            //direction_y = TMath::Sqrt(vy2*vy2+1*1);
            //double eta2_y = TMath::ATanH(1/direction_y); 

            double delta_eta_x = TMath::Abs(eta2_x-eta1_x); 
            double delta_eta_y = TMath::Abs(eta2_y-eta1_y); 
            if(false){
                if(delta_eta_x > eta_threshold(10*(clusterFirst->global().x() - collision_x), 0.010)){
                    hFilter->Fill(1);
                    continue;
                }
                if(delta_eta_y > eta_threshold(10*(clusterFirst->global().y() - collision_y), 0.010)){
                    hFilter->Fill(1);
                    continue;
                }
            }
            double vx = (clusterLast->global().x() - clusterFirst->global().x())/(clusterLast->global().z() - clusterFirst->global().z());
            double vy = (clusterLast->global().y() - clusterFirst->global().y())/(clusterLast->global().z() - clusterFirst->global().z());

            direction_tot = TMath::Sqrt(vx*vx+vy*vy+1*1);
            double eta = TMath::Log((direction_tot+1)/(direction_tot-1))/2.;
            direction_tot = TMath::Sqrt(vx*vx+1*1);
            double eta_x = TMath::Log((direction_tot+1)/(direction_tot-1))/2.;
            direction_tot = TMath::Sqrt(vy*vy+1*1);
            double eta_y = TMath::Log((direction_tot+1)/(direction_tot-1))/2.;
            if(apply_selections_){
                if(eta_x < min_eta_x_ || eta_y < min_eta_y_ ){
                    hFilter->Fill(2);
                    continue;
                }
                if(eta1_x < min_eta_x_ || eta1_y < min_eta_y_ ){
                    hFilter->Fill(2);
                    continue;
                }
                /*
                if(((clusterFirst->global().y() - collision_y)/(clusterLast->global().y() - clusterFirst->global().y()))<0){
                    continue;
                }
                if(((clusterFirst->global().x() - collision_x)/(clusterLast->global().x() - clusterFirst->global().x()))<0){
                    continue;
                }
                */
                if(eta > max_eta_){
                    hFilter->Fill(3);
                    continue;
                }
            }
            
            refTrack.addCluster(clusterFirst.get());
            refTrack.addCluster(clusterLast.get());
            refTrack.registerPlane(reference_first->getName(),
                                   reference_first->displacement().z(),
                                   reference_first->materialBudget(),
                                   reference_first->toLocal());
            refTrack.registerPlane(reference_last->getName(),
                                   reference_last->displacement().z(),
                                   reference_last->materialBudget(),
                                   reference_last->toLocal());

            refTrack.fit();

            // Now look for the spatially closest cluster on the next plane
            //auto interceptPointPV = refTrack.getIntercept(100);
            
            const std::string dummy_str2;
            auto state2 = refTrack.getIntercept(0);
            auto direction2 = refTrack.getDirection(dummy_str2);
            double interceptXPV = direction2.X()*100+state2.X();
            double interceptYPV = direction2.Y()*100+state2.Y();
            double distanceXPV = collision_x - interceptXPV;
            double distanceYPV = collision_y - interceptYPV- 1.5;

            //LOG(DEBUG) << "Cluster outside the cuts. Normalized distance: " << normPV;
            //LOG(DEBUG) << "Pb x:  " << collision_x;
            //LOG(DEBUG) << "Trk x: " << interceptXPV;
            //LOG(DEBUG) << "Pb x:  " << collision_y;
            //LOG(DEBUG) << "Trk y: " << interceptYPV;

            double normPV = (distanceXPV * distanceXPV) + (distanceYPV * distanceYPV);

            if(normPV > 0.5*0.5) {
                continue;
            }
            // Make a new track
            auto track = Track::Factory("straightline");
            track->setDistFromPV(normPV);
            track->addCluster(clusterFirst.get());
            track->addCluster(clusterLast.get());

            // Fit initial trajectory guess
            //LOG(DEBUG) << "Ref fitted";
            // Loop over each subsequent plane and look for a cluster within the timing cuts
            size_t detector_nr = 2;
            // Get all detectors here to also include passive layers which might contribute to scattering
            for(auto& detector : get_detectors()) {
                if(detector->isAuxiliary()) {
                    continue;
                }
                auto detectorID = detector->getName();
                LOG(TRACE) << "Registering detector " << detectorID << " at z = " << detector->displacement().z();

                // Add plane to track and trigger re-fit:
                refTrack.updatePlane(
                    detectorID, detector->displacement().z(), detector->materialBudget(), detector->toLocal());
                track->registerPlane(
                    detectorID, detector->displacement().z(), detector->materialBudget(), detector->toLocal());

                if(detector == reference_first || detector == reference_last) {
                    continue;
                }

                if(exclude_DUT_ && detector->isDUT()) {
                    //LOG(DEBUG) << "Skipping DUT plane.";
                    continue;
                }

                if(detector->isPassive()) {
                    //LOG(DEBUG) << "Skipping passive plane.";
                    continue;
                }

                // Determine whether a track can still be assembled given the number of current hits and the number of
                // detectors to come. Reduces computing time.
                detector_nr++;
                if(refTrack.getNClusters() + (trees.size() - detector_nr + 1) < min_hits_on_track_) {
                    //LOG(DEBUG) << "No chance to find a track - too few detectors left: " << refTrack.getNClusters() << " + " << trees.size() << " - " << detector_nr << " < " << min_hits_on_track_;
                    continue;
                }

                if(trees.count(detector) == 0) {
                    LOG(TRACE) << "Skipping detector " << detector->getName() << " as it has 0 clusters.";
                    continue;
                }

                // Get all neighbors within the timing cut
                //LOG(DEBUG) << "Searching for neighboring cluster on device " << detector->getName();
                ////LOG(DEBUG) << "- reference time is " << Units::display(refTrack.timestamp(), {"ns", "us", "s"});
                Cluster* closestCluster = nullptr;

                // Use spatial cut only as initial value (check if cluster is ellipse defined by cuts is done below):
                double closestClusterDistance = sqrt(spatial_cuts_[detectorID] * spatial_cuts_[detectorID] +
                                                     spatial_cuts_[detectorID] * spatial_cuts_[detectorID]);

                auto neighbors = trees[detector].getAllElementsInTimeWindow(refTrack.timestamp(), 1e99);

                //LOG(DEBUG) << "- found " << neighbors.size() << " neighbors within the correct time window on " << detectorID;

                // Now look for the spatially closest cluster on the next plane
                PositionVector3D<Cartesian3D<double>> interceptPoint = detector->getLocalIntercept(&refTrack);
                double interceptX = interceptPoint.X();
                double interceptY = interceptPoint.Y();

                for(size_t ne = 0; ne < neighbors.size(); ne++) {
                    auto newCluster = neighbors[ne].get();

                    // Calculate the distance to the previous plane's cluster/intercept
                    double distanceX = interceptX - newCluster->local().x();
                    double distanceY = interceptY - newCluster->local().y();
                    double distance = sqrt(distanceX * distanceX + distanceY * distanceY);

                    // Check if newCluster lies within ellipse defined by spatial cuts around intercept,
                    // following this example:
                    // https://www.geeksforgeeks.org/check-if-a-point-is-inside-outside-or-on-the-ellipse/
                    //
                    // ellipse defined by: x^2/a^2 + y^2/b^2 = 1: on ellipse,
                    //                                       > 1: outside,
                    //                                       < 1: inside
                    // Continue if outside of ellipse:

                    double norm = (distanceX * distanceX) / (spatial_cuts_[detectorID] * spatial_cuts_[detectorID]) +
                                  (distanceY * distanceY) / (spatial_cuts_[detectorID] * spatial_cuts_[detectorID]);

                    if(norm > 1) {
                        //LOG(DEBUG) << "Cluster outside the cuts. Normalized distance: " << norm;
                        continue;
                    }

                    // If this is the closest keep it for now
                    if(distance < closestClusterDistance) {
                        closestClusterDistance = distance;
                        closestCluster = newCluster;
                    }
                }

                if(closestCluster == nullptr) {
                    //LOG(DEBUG) << "No cluster within spatial cut";
                    continue;
                }

                // Add the cluster to the track
                refTrack.addCluster(closestCluster);
                track->addCluster(closestCluster);

                ////LOG(DEBUG) << "- added cluster to track";
            }

            // check if track has required detector(s):
            auto foundRequiredDetector = [this](Track* t) {
                for(auto& requireDet : require_detectors_) {
                    if(!requireDet.empty() && !t->hasDetector(requireDet)) {
                        //LOG(DEBUG) << "No cluster from required detector " << requireDet << " on the track.";
                        return false;
                    }
                }
                return true;
            };
            if(!foundRequiredDetector(track.get())) {
                continue;
            }

            // Now should have a track with one cluster from each plane
            if(track->getNClusters() < min_hits_on_track_) {
                //LOG(DEBUG) << "Not enough clusters on the track, found " << track->getNClusters() << " but " << min_hits_on_track_ << " required.";
                continue;
            }

            // Fit the track
            track->fit();

            if(reject_by_ROI_ && track->isFitted()) {
                // check if the track is within ROI for all detectors
                auto ds = get_regular_detectors(!exclude_DUT_);
                auto out_of_roi =
                    std::find_if(ds.begin(), ds.end(), [track](const auto& d) { return !d->isWithinROI(track.get()); });
                if(out_of_roi != ds.end()) {
                    //LOG(DEBUG) << "Rejecting track outside of ROI of detector " << out_of_roi->get()->getName();
                    continue;
                }
            }
            // save the track
            if(track->isFitted()) {
                tracks.push_back(track);
            } else {
                LOG_N(WARNING, 100) << "Rejected a track due to failure in fitting";
                continue;
            }
        }
    }

    auto duplicated_hit = [this](const Track* a, const Track* b) {
        for(auto d : get_regular_detectors(!exclude_DUT_)) {
            if(a->getClusterFromDetector(d->getName()) == b->getClusterFromDetector(d->getName()) &&
               !(b->getClusterFromDetector(d->getName()) == nullptr)) {
                ////LOG(DEBUG) << "Duplicated hit on " << d->getName() << ": rejecting track";
                return true;
            }
        }
        return false;
    };
    // Save the tracks on the clipboard
    if(tracks.size() > 0) {

        // if requested ensure unique usage of clusters
        if(unique_cluster_usage_ && tracks.size() > 1) {
            // sort by chi2:
            LOG_ONCE(WARNING) << "Rejecting tracks with same hits";
            /*
            std::sort(tracks.begin(), tracks.end(), [](const shared_ptr<Track> a, const shared_ptr<Track> b) {
                return (a->getChi2() / static_cast<double>(a->getNdof())) <
                       (b->getChi2() / static_cast<double>(b->getNdof()));
            });
            */            
            std::sort(tracks.begin(), tracks.end(), [](const shared_ptr<Track> a, const shared_ptr<Track> b) {
                return a->getDistFromPV() < b->getDistFromPV();
            });
            // remove tracks with hit that is used twice
            auto track1 = tracks.begin();
            while(track1 != tracks.end()) {
                auto track2 = track1 + 1;
                while(track2 != tracks.end()) {
                    // if hit is used twice delete the track
                    if(duplicated_hit(track2->get(), track1->get())) {
                        track2 = tracks.erase(track2);
                    } else {
                        track2++;
                    }
                }
                track1++;
            }
        }
        clipboard->putData(tracks);
    }
    for(auto track : tracks) {
        // Fill track time within event (relative to event start)
        auto event = clipboard->getEvent();

        trackChi2->Fill(track->getChi2());
        clustersPerTrack->Fill(static_cast<double>(track->getNClusters()));
        trackChi2ndof->Fill(track->getChi2ndof());
        trackAngleX->Fill(atan(track->getDirection(track->getClusters().front()->detectorID()).X()));
        trackAngleY->Fill(atan(track->getDirection(track->getClusters().front()->detectorID()).Y()));
        
        const std::string dummy_str;
        auto state = track->getIntercept(0);
        auto direction = track->getDirection(dummy_str);

        CorryTrack trackCorry;
        trackCorry.vx = static_cast<float>(direction.X());
        trackCorry.vy = static_cast<float>(direction.Y());
        trackCorry.px = static_cast<float>(state.X());
        trackCorry.py = static_cast<float>(state.Y());
        trackCorry.chi2ndof = static_cast<float>(track->getChi2ndof());

        trackCorry.size1 = 0;
        trackCorry.size2 = 0;
        trackCorry.size3 = 0;
        trackCorry.size4 = 0;
        trackCorry.size5 = 0;
        trackCorry.size6 = 0;

        // Make residuals
        auto trackClusters = track->getClusters();
        for(auto& trackCluster : trackClusters) {
            string detectorID = trackCluster->detectorID();
            
            string str1 ("2");
            string str2 ("3");
            string str3 ("4");
            string str4 ("5");
            string str5 ("6");
            if (detectorID.find(str1) != string::npos)
                trackCorry.size1 = static_cast<int>(trackCluster->size());
            else if (detectorID.find(str2) != string::npos)
                trackCorry.size2 = static_cast<int>(trackCluster->size());
            else if (detectorID.find(str3) != string::npos)
                trackCorry.size3 = static_cast<int>(trackCluster->size());
            else if (detectorID.find(str4) != string::npos)
                trackCorry.size4 = static_cast<int>(trackCluster->size());
            else if (detectorID.find(str5) != string::npos)
                trackCorry.size5 = static_cast<int>(trackCluster->size());
            else
                trackCorry.size6 = static_cast<int>(trackCluster->size());
            
            ROOT::Math::XYZPoint globalRes = track->getGlobalResidual(detectorID);
            ROOT::Math::XYPoint localRes = track->getLocalResidual(detectorID);

            residualsX_local[detectorID]->Fill(localRes.X());
            residualsX_global[detectorID]->Fill(globalRes.X());
            residualsX_vs_positionX_global[detectorID]->Fill(globalRes.X(), trackCluster->global().x());
            residualsX_vs_positionY_global[detectorID]->Fill(globalRes.X(), trackCluster->global().y());

            pullX_local[detectorID]->Fill(localRes.x() / track->getClusterFromDetector(detectorID)->errorX());
            pullX_global[detectorID]->Fill(globalRes.x() / track->getClusterFromDetector(detectorID)->errorX());

            pullY_local[detectorID]->Fill(localRes.Y() / track->getClusterFromDetector(detectorID)->errorY());
            pullY_global[detectorID]->Fill(globalRes.Y() / track->getClusterFromDetector(detectorID)->errorY());

            if(trackCluster->columnWidth() == 1) {
                residualsXwidth1_local[detectorID]->Fill(localRes.X());
                residualsXwidth1_global[detectorID]->Fill(globalRes.X());
            } else if(trackCluster->columnWidth() == 2) {
                residualsXwidth2_local[detectorID]->Fill(localRes.X());
                residualsXwidth2_global[detectorID]->Fill(globalRes.X());
            } else if(trackCluster->columnWidth() == 3) {
                residualsXwidth3_local[detectorID]->Fill(localRes.X());
                residualsXwidth3_global[detectorID]->Fill(globalRes.X());
            }

            residualsY_local[detectorID]->Fill(localRes.Y());
            residualsY_global[detectorID]->Fill(globalRes.Y());
            residualsY_vs_positionY_global[detectorID]->Fill(globalRes.Y(), trackCluster->global().y());
            residualsY_vs_positionX_global[detectorID]->Fill(globalRes.Y(), trackCluster->global().x());

            if(trackCluster->rowWidth() == 1) {
                residualsYwidth1_local[detectorID]->Fill(localRes.Y());
                residualsYwidth1_global[detectorID]->Fill(globalRes.Y());
            } else if(trackCluster->rowWidth() == 2) {
                residualsYwidth2_local[detectorID]->Fill(localRes.Y());
                residualsYwidth2_global[detectorID]->Fill(globalRes.Y());
            } else if(trackCluster->rowWidth() == 3) {
                residualsYwidth3_local[detectorID]->Fill(localRes.Y());
                residualsYwidth3_global[detectorID]->Fill(globalRes.Y());
            }
        }

        for(auto& detector : get_regular_detectors(true)) {
            auto det = detector->getName();

            auto local = detector->getLocalIntercept(track.get());
            auto row = detector->getRow(local);
            auto col = detector->getColumn(local);
            LOG(TRACE) << "Local col/row intersect of track: " << col << "\t" << row;
            local_intersects_[det]->Fill(col, row);

            auto global = detector->getIntercept(track.get());
            global_intersects_[det]->Fill(global.X(), global.Y());
        }

        event_corry->tracks.push_back(trackCorry);
    }
    tracksPerEvent->Fill(static_cast<double>(tracks.size()));

    //LOG(DEBUG) << "Fill Tree";
    event_info->Fill();
    //LOG(DEBUG) << "clear";
    event_corry->tracks.clear();

    //LOG(DEBUG) << "End of event";
    return StatusCode::Success;
}
