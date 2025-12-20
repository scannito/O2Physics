// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \file PhiMesonCandidateProducer.cxx
/// \brief Table producer for  Phi mesons candidates
/// \author Stefano Cannito (stefano.cannito@cern.ch)

#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Tasks/Strangeness/testTable.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/Core/TableHelper.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <Math/Vector4D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THn.h>
#include <TMath.h>
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct PhiMesonCandProducer {
  // Produce the table with the event selection information
  Produces<aod::PhimesonCandidate> phimesonCandidate;

  HistogramRegistry histos{"histoQA", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurables for phi's daughter tracks selection
  struct : ConfigurableGroup {
    Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0f, "Cut on charge"};
    Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"};
    Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};
    Configurable<float> cMinKaonPtcut{"cMinKaonPtcut", 0.15f, "Track minimum pt cut"};
    Configurable<float> etaMax{"etaMax", 0.8f, "eta max"};
    Configurable<float> pTToUseTOF{"pTToUseTOF", 0.5f, "pT above which use TOF"};
    Configurable<float> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0f, "Track DCAz cut to PV Maximum"};
    Configurable<std::vector<float>> cMaxDCArToPVPhi{"cMaxDCArToPVPhi", {0.004f, 0.013f, 1.0f}, "Track DCAr cut to PV for Phi"};

    Configurable<bool> cfgIsDCAzParameterized{"cfgIsDCAzParameterized", false, "IsDCAzParameterized"};

    Configurable<float> nSigmaCutTPCKa{"nSigmaCutTPCKa", 3.0f, "Value of the TPC Nsigma cut for Kaons"};
    Configurable<float> nSigmaCutCombinedKa{"nSigmaCutCombinedKa", 3.0f, "Value of the TOF Nsigma cut for Kaons"};

    Configurable<int> minTPCnClsFound{"minTPCnClsFound", 70, "min number of found TPC clusters"};
    Configurable<int> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70, "min number of TPC crossed rows"};
    Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
    Configurable<int> minITSnCls{"minITSnCls", 4, "min number of ITS clusters"};
    Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
  } trackConfigs;

  // Configurables on phi selection
  struct : ConfigurableGroup {
    Configurable<float> lowMPhi{"lowMPhi", 1.0095f, "Upper limits on Phi mass for signal extraction"};
    Configurable<float> upMPhi{"upMPhi", 1.029f, "Upper limits on Phi mass for signal extraction"};

    Configurable<float> minPhiPt{"minPhiPt", 0.4f, "Minimum pT for Phi"};
    Configurable<float> maxPhiPt{"maxPhiPt", 10.0f, "Maximum pT for Phi"};

    Configurable<float> cfgYAcceptance{"cfgYAcceptance", 0.5f, "Rapidity acceptance"};
  } phiConfigs;

  // Defining the type of the collisions for data and MC
  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>;
  using SimCollisions = soa::Join<SelCollisions, aod::McCollisionLabels>;
  using MCCollisions = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;

  // Defining the type of the phi's daughter tracks for data and MC
  using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>;
  using FullMCTracks = soa::Join<FullTracks, aod::McTrackLabels>;

  Partition<FullTracks> posTracks = aod::track::signed1Pt > trackConfigs.cfgCutCharge;
  Partition<FullTracks> negTracks = aod::track::signed1Pt < trackConfigs.cfgCutCharge;

  Partition<FullMCTracks> posMCTracks = aod::track::signed1Pt > trackConfigs.cfgCutCharge;
  Partition<FullMCTracks> negMCTracks = aod::track::signed1Pt < trackConfigs.cfgCutCharge;

  // Manual slicing
  Preslice<aod::Tracks> trackPerCollision = aod::track::collisionId;
  SliceCache cache;

  // Constants
  double massKa = o2::constants::physics::MassKPlus;

  void init(InitContext&)
  {
    // Defining histogram axes
    AxisSpec vertexZAxis = {100, -cutZVertex, cutZVertex, "vrtx_{Z} [cm]"};
    AxisSpec multAxis = {120, 0.0f, 120.0f, "centFT0M"};
    AxisSpec binnedmultAxis{(std::vector<double>)binsMult, "centFT0M"};

    // Booking histograms for event selection QA
    // Number of events per selection in Data
    histos.add("hEventSelectionData", "hEventSelectionData", kTH1F, {{5, -0.5f, 4.5f}});
    histos.get<TH1>(HIST("hEventSelectionData"))->GetXaxis()->SetBinLabel(1, "All collisions");
    histos.get<TH1>(HIST("hEventSelectionData"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
    histos.get<TH1>(HIST("hEventSelectionData"))->GetXaxis()->SetBinLabel(3, "posZ cut");
    histos.get<TH1>(HIST("hEventSelectionData"))->GetXaxis()->SetBinLabel(4, "INEL>0 cut");
    histos.get<TH1>(HIST("hEventSelectionData"))->GetXaxis()->SetBinLabel(5, "With at least a #phi cand");

    // Number of MC events per selection in MC
    histos.add("hEventSelectionMC", "hEventSelectionMC", kTH1F, {{8, -0.5f, 7.5f}});
    histos.get<TH1>(HIST("hEventSelectionMC"))->GetXaxis()->SetBinLabel(1, "All collisions");
    histos.get<TH1>(HIST("hEventSelectionMC"))->GetXaxis()->SetBinLabel(2, "kIsTriggerTVX");
    histos.get<TH1>(HIST("hEventSelectionMC"))->GetXaxis()->SetBinLabel(3, "kNoTimeFrameBorder");
    histos.get<TH1>(HIST("hEventSelectionMC"))->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
    histos.get<TH1>(HIST("hEventSelectionMC"))->GetXaxis()->SetBinLabel(5, "posZ cut");
    histos.get<TH1>(HIST("hEventSelectionMC"))->GetXaxis()->SetBinLabel(6, "INEL>0 cut");
    histos.get<TH1>(HIST("hEventSelectionMC"))->GetXaxis()->SetBinLabel(7, "With at least a gen coll");
    histos.get<TH1>(HIST("hEventSelectionMC"))->GetXaxis()->SetBinLabel(8, "With at least a #phi cand");

    // Event information
    histos.add("hVertexZ", "Vertex Z", kTH1F, {vertexZAxis});
    histos.add("hVertexZWPhi", "Vertex Z with a Phi Candidate", kTH1F, {vertexZAxis});
    histos.add("hMultiplicityPercent", "Multiplicity Percentile", kTH1F, {multAxis});
    histos.add("hMultiplicityPercentWPhi", "Multiplicity Percentile in Events with a Phi Candidate", kTH1F, {multAxis});
    histos.add("h2VertexZvsMult", "Vertex Z vs Multiplicity Percentile", kTH2F, {vertexZAxis, binnedmultAxis});
    histos.add("h2VertexZvsMultWPhi", "Vertex Z vs Multiplicity Percentile with a Phi Candidate", kTH2F, {vertexZAxis, binnedmultAxis});

    // Phi's daughter tracks information
    /*histos.add("hEta", "Eta of Kaon candidates", kTH1F, {{100, -1.0f, 1.0f}});
    histos.add("hNsigmaKaonTPC", "NsigmaKaon TPC distribution vs pt", kTH2F, {{100, 0.0, 5.0, "#it{p} (GeV/#it{c})"}, {100, -10.0f, 10.0f}});
    histos.add("hNsigmaKaonTOF", "NsigmaKaon TOF distribution vs pt", kTH2F, {{100, 0.0, 5.0, "#it{p} (GeV/#it{c})"}, {100, -10.0f, 10.0f}});
    histos.add("h2DauTracksPhiDCAxy", "DCAxy distribution vs pt", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    histos.add("h2DauTracksPhiDCAz", "DCAz distribution vs pt", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});*/
  }

  // Event selection and QA filling
  template <bool isMC, typename MC_T = void, typename T>
  bool defaultEventSelection(const T& collision)
  {
    float multPercentile = 0.0f;

    if constexpr (!isMC) {                         // data event
      histos.fill(HIST("hEventSelectionData"), 0); // all collisions
      if (!collision.sel8())
        return false;
      histos.fill(HIST("hEventSelectionData"), 1); // sel8 collisions
      if (std::abs(collision.posZ()) >= cutZVertex)
        return false;
      histos.fill(HIST("hEventSelectionData"), 2); // vertex-Z selected
      if (!collision.isInelGt0())
        return false;
      histos.fill(HIST("hEventSelectionData"), 3); // INEL>0 collisions

      multPercentile = collision.centFT0M();
    } else { // MCreco event
      static_assert(!std::is_same_v<MC_T, void>, "Need to set MC_T to MCCollisions for isMC = true");

      histos.fill(HIST("hEventSelectionMC"), 0); // all collisions
      if (!collision.selection_bit(aod::evsel::kIsTriggerTVX))
        return false;
      histos.fill(HIST("hEventSelectionMC"), 1); // kIsTriggerTVX collisions
      if (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
        return false;
      histos.fill(HIST("hEventSelectionMC"), 2); // kNoTimeFrameBorder collisions
      if (cfgiskNoITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))
        return false;
      histos.fill(HIST("hEventSelectionMC"), 3); // kNoITSROFrameBorder collisions (by default not requested by the selection)
      if (std::abs(collision.posZ()) > cutZVertex)
        return false;
      histos.fill(HIST("hEventSelectionMC"), 4); // vertex-Z selected
      if (!collision.isInelGt0())
        return false;
      histos.fill(HIST("hEventSelectionMC"), 5); // INEL>0 collisions
      if (!collision.has_mcCollision())
        return false;
      histos.fill(HIST("hEventSelectionMC"), 6); // with at least a gen collision

      const auto& mcCollision = collision.template mcCollision_as<MC_T>();
      multPercentile = mcCollision.centFT0M();
    }

    histos.fill(HIST("hVertexZ"), collision.posZ());
    histos.fill(HIST("hMultiplicityPercent"), multPercentile);
    histos.fill(HIST("h2VertexZvsMult"), collision.posZ(), multPercentile);

    return true;
  }

  // Topological track selection
  template <typename T>
  bool selectionTrackResonance(const T& track)
  {
    if (trackConfigs.cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (trackConfigs.cfgPVContributor && !track.isPVContributor())
      return false;

    if (track.tpcNClsFound() < trackConfigs.minTPCnClsFound)
      return false;

    if (track.pt() < trackConfigs.cMinKaonPtcut)
      return false;

    if (std::abs(track.dcaXY()) > trackConfigs.cMaxDCArToPVPhi->at(0) + (trackConfigs.cMaxDCArToPVPhi->at(1) / std::pow(track.pt(), trackConfigs.cMaxDCArToPVPhi->at(2))))
      return false;
    if (std::abs(track.dcaZ()) > trackConfigs.cMaxDCAzToPVcut)
      return false;

    return true;
  }

  // PIDKaon track selection
  template <typename T>
  bool selectionPIDKaonpTdependent(const T& track)
  {
    if (track.pt() < trackConfigs.pTToUseTOF && std::abs(track.tpcNSigmaKa()) < trackConfigs.nSigmaCutTPCKa)
      return true;
    if (track.pt() >= trackConfigs.pTToUseTOF && track.hasTOF() && (std::pow(track.tofNSigmaKa(), 2) + std::pow(track.tpcNSigmaKa(), 2)) < std::pow(trackConfigs.nSigmaCutCombinedKa, 2))
      return true;

    return false;
  }

  // Reconstruct the Phi candidate
  template <typename T1, typename T2>
  ROOT::Math::PxPyPzMVector recMother(const T1& track1, const T2& track2, float masscand1, float masscand2)
  {
    ROOT::Math::PxPyPzMVector daughter1(track1.px(), track1.py(), track1.pz(), masscand1); // set the daughter1 4-momentum
    ROOT::Math::PxPyPzMVector daughter2(track2.px(), track2.py(), track2.pz(), masscand2); // set the daughter2 4-momentum
    ROOT::Math::PxPyPzMVector mother = daughter1 + daughter2;                              // calculate the mother 4-momentum

    return mother;
  }

  template <bool isMC, typename MC_T = void, typename T>
  bool eventHasPhi(const T& collision)
  {
    auto countPhiCandidates = [&](auto& pTracks, auto& nTracks) {
      // auto posThisColl = pTracks.sliceBy(trackPerCollision, collision.globalIndex());
      // auto negThisColl = nTracks.sliceBy(trackPerCollision, collision.globalIndex());
      auto posThisColl = pTracks.sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      auto negThisColl = nTracks.sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

      int count = 0;

      for (const auto& track1 : posThisColl) {
        if (!selectionTrackResonance(track1) || !selectionPIDKaonpTdependent(track1))
          continue;

        for (const auto& track2 : negThisColl) {
          if (!selectionTrackResonance(track2) || !selectionPIDKaonpTdependent(track2))
            continue;

          // histos.fill(HIST("hEta"), track1.eta());
          // histos.fill(HIST("hNsigmaKaonTPC"), track1.tpcInnerParam(), track1.tpcNSigmaKa());
          // histos.fill(HIST("hNsigmaKaonTOF"), track1.tpcInnerParam(), track1.tofNSigmaKa());
          // histos.fill(HIST("h2DauTracksPhiDCAxy"), track1.pt(), track1.dcaXY());
          // histos.fill(HIST("h2DauTracksPhiDCAz"), track1.pt(), track1.dcaZ());

          ROOT::Math::PxPyPzMVector recPhi = recMother(track1, track2, massKa, massKa);

          if (recPhi.Pt() < phiConfigs.minPhiPt)
            continue;
          if (recPhi.M() < phiConfigs.lowMPhi || recPhi.M() > phiConfigs.upMPhi)
            continue;
          if (std::abs(recPhi.Rapidity()) > phiConfigs.cfgYAcceptance)
            continue;

          count++;
        }
      }

      return count;
    };

    float multPercentile = 0.0f;

    if constexpr (!isMC) {
      if (countPhiCandidates(posTracks, negTracks) == 0)
        return false;
      histos.fill(HIST("hEventSelectionData"), 4);

      multPercentile = collision.centFT0M();
    } else {
      static_assert(!std::is_same_v<MC_T, void>, "Need to set MC_T to MCCollisions for isMC = true");

      if (countPhiCandidates(posMCTracks, negMCTracks) == 0)
        return false;
      histos.fill(HIST("hEventSelectionMC"), 7);

      const auto& mcCollision = collision.template mcCollision_as<MC_T>();
      multPercentile = mcCollision.centFT0M();
    }

    histos.fill(HIST("hVertexZWPhi"), collision.posZ());
    histos.fill(HIST("hMultiplicityPercentWPhi"), multPercentile);
    histos.fill(HIST("h2VertexZvsMultWPhi"), collision.posZ(), multPercentile);

    return true;
  }

  void processData(SelCollisions::iterator const& collision, FullTracks const&)
  {
    phimesonSelection(defaultEventSelection<false>(collision) && eventHasPhi<false>(collision));
  }

  PROCESS_SWITCH(PhiMesonSelCollision, processData, "Process function to select events with Phi mesons in Data", true);

  void processMC(SimCollisions::iterator const& collision, MCCollisions const&, FullMCTracks const&)
  {
    phimesonSelection(defaultEventSelection<true, MCCollisions>(collision) && eventHasPhi<true, MCCollisions>(collision));
  }

  PROCESS_SWITCH(PhiMesonSelCollision, processMC, "Process function to select events with Phi mesons in MC", false);

  /*
  template <bool isMC, typename MC_T = void, typename T1, typename T2>
  bool eventHasPhi(const T1& collision, const T2& posTracks, const T2& negTracks)
  {
    float multPercentile = 0.0f;

    int nPhi = 0;

    for (const auto& track1 : posTracks) {
      if (!selectionTrackResonance(track1) || !selectionPIDKaonpTdependent(track1))
        continue;

      for (const auto& track2 : negTracks) {
        if (!selectionTrackResonance(track2) || !selectionPIDKaonpTdependent(track2))
          continue;

        // histos.fill(HIST("hEta"), track1.eta());
        // histos.fill(HIST("hNsigmaKaonTPC"), track1.tpcInnerParam(), track1.tpcNSigmaKa());
        // histos.fill(HIST("hNsigmaKaonTOF"), track1.tpcInnerParam(), track1.tofNSigmaKa());
        // histos.fill(HIST("h2DauTracksPhiDCAxy"), track1.pt(), track1.dcaXY());
        // histos.fill(HIST("h2DauTracksPhiDCAz"), track1.pt(), track1.dcaZ());

        ROOT::Math::PxPyPzMVector recPhi = recMother(track1, track2, massKa, massKa);

        if (recPhi.Pt() < phiConfigs.minPhiPt)
          continue;
        if (recPhi.M() < phiConfigs.lowMPhi || recPhi.M() > phiConfigs.upMPhi)
          continue;
        if (std::abs(recPhi.Rapidity()) > phiConfigs.cfgYAcceptance)
          continue;

        nPhi++;
      }
    }

    if (nPhi == 0)
      return false;

    if constexpr (!isMC) {
      histos.fill(HIST("hEventSelectionData"), 4);

      multPercentile = collision.centFT0M();
    } else {
      static_assert(!std::is_same_v<MC_T, void>, "Need to set MC_T to MCCollisions for isMC = true");

      histos.fill(HIST("hEventSelectionMC"), 7);

      const auto& mcCollision = collision.template mcCollision_as<MC_T>();
      multPercentile = mcCollision.centFT0M();
    }

    histos.fill(HIST("hVertexZWPhi"), collision.posZ());
    histos.fill(HIST("hMultiplicityPercentWPhi"), multPercentile);
    histos.fill(HIST("h2VertexZvsMultWPhi"), collision.posZ(), multPercentile);

    return true;
  }

  void processData(SelCollisions::iterator const& collision, FullTracks const&)
  {
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    phimesonSelection(defaultEventSelection<false>(collision) && eventHasPhi<false>(collision, posThisColl, negThisColl));
  }

  PROCESS_SWITCH(PhiMesonSelCollision, processData, "Process function to select events with Phi mesons in Data", true);

  void processMC(SimCollisions::iterator const& collision, MCCollisions const&, FullMCTracks const&)
  {
    auto posThisColl = posMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    phimesonSelection(defaultEventSelection<true, MCCollisions>(collision) && eventHasPhi<true, MCCollisions>(collision, posThisColl, negThisColl));
  }

  PROCESS_SWITCH(PhiMesonSelCollision, processMC, "Process function to select events with Phi mesons in MC", false);
  */
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PhiMesonSelCollision>(cfgc)};
}
