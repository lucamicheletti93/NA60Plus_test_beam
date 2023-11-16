/**
 * @file
 * @brief Definition of module Tracking4D
 *
 * @copyright Copyright (c) 2017-2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef TRACKING4D_H
#define TRACKING4D_H 1

#include <TCanvas.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <iostream>
#include <vector>

#include "core/module/Module.hpp"
#include "objects/Cluster.hpp"
#include "objects/Pixel.hpp"
#include "objects/Track.hpp"
#include "objects/CorryTree.hpp"

namespace corryvreckan {
    /** @ingroup Modules
     */
    class Tracking4D : public Module {

    public:
        // Constructors and destructors
        Tracking4D(Configuration& config, std::vector<std::shared_ptr<Detector>> detectors);
        ~Tracking4D() {}

        // Functions
        void initialize() override;
        StatusCode run(const std::shared_ptr<Clipboard>& clipboard) override;

    private:
        // Histograms
        TH1F* trackChi2;
        TH1F* clustersPerTrack;
        TH1F* trackChi2ndof;
        TH1F* tracksPerEvent;
        TH1F* trackAngleX;
        TH1F* trackAngleY;
        TH1F* hFilter;
        std::map<std::string, TH1F*> residualsX_local;
        std::map<std::string, TH1F*> residualsXwidth1_local;
        std::map<std::string, TH1F*> residualsXwidth2_local;
        std::map<std::string, TH1F*> residualsXwidth3_local;
        std::map<std::string, TH1F*> pullY_local;
        std::map<std::string, TH1F*> residualsY_local;
        std::map<std::string, TH1F*> residualsYwidth1_local;
        std::map<std::string, TH1F*> residualsYwidth2_local;
        std::map<std::string, TH1F*> residualsYwidth3_local;
        std::map<std::string, TH1F*> pullX_local;

        std::map<std::string, TH1F*> residualsX_global;
        std::map<std::string, TH2F*> residualsX_vs_positionX_global;
        std::map<std::string, TH2F*> residualsX_vs_positionY_global;
        std::map<std::string, TH1F*> residualsXwidth1_global;
        std::map<std::string, TH1F*> residualsXwidth2_global;
        std::map<std::string, TH1F*> residualsXwidth3_global;
        std::map<std::string, TH1F*> pullX_global;
        std::map<std::string, TH1F*> residualsY_global;
        std::map<std::string, TH2F*> residualsY_vs_positionY_global;
        std::map<std::string, TH2F*> residualsY_vs_positionX_global;
        std::map<std::string, TH1F*> residualsYwidth1_global;
        std::map<std::string, TH1F*> residualsYwidth2_global;
        std::map<std::string, TH1F*> residualsYwidth3_global;
        std::map<std::string, TH1F*> pullY_global;

        std::map<std::string, TH2F*> local_intersects_;
        std::map<std::string, TH2F*> global_intersects_;

        // Cuts for tracking
        double max_plot_chi2_;
        size_t min_hits_on_track_;
        bool exclude_DUT_;
        bool apply_selections_;
        bool reject_by_ROI_;
        bool unique_cluster_usage_;
        float max_delta_eta_;
        float min_eta_y_;
        float min_eta_x_;
        float max_eta_;
        float target_position_;

        std::map<std::string, float> spatial_cuts_;
        float spatial_cuts_max_;
        int max_range_multiplicity_;
        std::vector<std::string> require_detectors_;
        std::vector<std::string> exclude_from_seed_;
        std::map<std::shared_ptr<Detector>, XYVector> spatial_cuts_tmp_;



        // Tree with the tracks and event information
        TTree* event_info;
        CorryEvent* event_corry;

        
        // Define your custom function
        double eta_threshold(double x, double params);
    };
} // namespace corryvreckan
#endif // TRACKING4D_H
