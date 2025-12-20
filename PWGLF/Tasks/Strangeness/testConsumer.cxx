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
/// \file PhiStrangenessCorrelation.cxx
/// \brief Analysis task for phi-strangeness rapidity correlations analysis
/// \author Stefano Cannito (stefano.cannito@cern.ch)

#include "PWGLF/DataModel/LFStrangenessTables.h"
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
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THn.h>
#include <TList.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TRandom.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct BoundEfficiencyMap {
  using CoordsTuple = std::tuple<float, float, float>;

  const TH3* effMap;
  CoordsTuple coords;

  BoundEfficiencyMap(const std::shared_ptr<TH3>& effMap, float x, float y, float z) : effMap(effMap.get()), coords(x, y, z) {}
  BoundEfficiencyMap(const std::shared_ptr<TH3>& effMap, const CoordsTuple& coords) : effMap(effMap.get()), coords(coords) {}

  float getEfficiency() const
  {
    if (!effMap) {
      return 1.0f;
    }

    const auto& [x, y, z] = coords;
    return effMap->Interpolate(x, y, z);
  }
};

struct PhiStrangenessCorrelation {
  HistogramRegistry histos{"phiStrangenessCorrelation", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable on multiplicity bins
  Configurable<std::vector<double>> binsMult{"binsMult", {0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0}, "Multiplicity bin limits"};

  // Configurables for tracks selection
  struct : ConfigurableGroup {
    Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0f, "Cut on charge"};
    Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"};
    Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};
    Configurable<float> cMinKaonPtcut{"cMinKaonPtcut", 0.15f, "Track minimum pt cut"};
    Configurable<float> etaMax{"etaMax", 0.8f, "eta max"};
    Configurable<float> pTToUseTOF{"pTToUseTOF", 0.5f, "pT above which use TOF"};
    Configurable<float> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0f, "Track DCAz cut to PV Maximum"};
    Configurable<std::vector<float>> cMaxDCArToPVPhi{"cMaxDCArToPVPhi", {0.004f, 0.013f, 1.0f}, "Track DCAr cut to PV for Phi"};

    Configurable<bool> cfgIsTOFChecked{"cfgIsTOFChecked", true, "Is TOF checked in PID for pions"};
    Configurable<std::vector<float>> cMaxDCArToPVPion{"cMaxDCArToPVPion", {0.004f, 0.013f, 1.0f}, "Track DCAr cut to PV for Pions"};
    Configurable<bool> cfgIsDCAzParameterized{"cfgIsDCAzParameterized", false, "IsDCAzParameterized"};
    Configurable<std::vector<float>> cMaxDCAzToPVPion{"cMaxDCAzToPVPion", {0.004f, 0.013f, 1.0f}, "Track DCAz cut to PV for Pions"};

    Configurable<float> nSigmaCutTPCKa{"nSigmaCutTPCKa", 3.0f, "Value of the TPC Nsigma cut for Kaons"};
    Configurable<float> nSigmaCutCombinedKa{"nSigmaCutCombinedKa", 3.0f, "Value of the TOF Nsigma cut for Kaons"};

    Configurable<float> nSigmaCutTPCPion{"nSigmaCutTPCPion", 4.0f, "Value of the TPC Nsigma cut for Pions"};
    Configurable<float> cMinPionPtcut{"cMinPionPtcut", 0.2f, "Track minimum pt cut"};

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
  } phiConfigs;

  // Configurables on phi pT bins
  Configurable<std::vector<double>> binspTPhi{"binspTPhi", {0.4, 0.8, 1.4, 2.0, 2.8, 4.0, 6.0, 10.0}, "pT bin limits for Phi"};

  // Configurables for V0 selection
  struct : ConfigurableGroup {
    Configurable<float> v0SettingCosPA{"v0SettingCosPA", 0.98f, "V0 CosPA"};
    Configurable<float> v0SettingRadius{"v0SettingRadius", 0.5f, "v0radius"};
    Configurable<float> v0SettingDCAV0Dau{"v0SettingDCAV0Dau", 1.0f, "DCA V0 Daughters"};
    Configurable<float> v0SettingDCAPosToPV{"v0SettingDCAPosToPV", 0.1f, "DCA Pos To PV"};
    Configurable<float> v0SettingDCANegToPV{"v0SettingDCANegToPV", 0.1f, "DCA Neg To PV"};
    Configurable<float> v0SettingMinPt{"v0SettingMinPt", 0.1f, "V0 min pt"};

    Configurable<bool> cfgFurtherV0Selection{"cfgFurtherV0Selection", false, "Further V0 selection"};
    Configurable<float> ctauK0s{"ctauK0s", 20.0f, "C tau K0s(cm)"};
    Configurable<float> paramArmenterosCut{"paramArmenterosCut", 0.2f, "parameter Armenteros Cut"};
    Configurable<float> v0rejK0s{"v0rejK0s", 0.005f, "V0 rej K0s"};

    Configurable<float> lowMK0S{"lowMK0S", 0.48f, "Lower limit on K0Short mass"};
    Configurable<float> upMK0S{"upMK0S", 0.52f, "Upper limit on K0Short mass"};
  } v0Configs;

  // Configurable on K0S pT bins
  Configurable<std::vector<double>> binspTK0S{"binspTK0S", {0.1, 0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0}, "pT bin limits for K0S"};

  // Configurable on pion pT bins
  Configurable<std::vector<double>> binspTPi{"binspTPi", {0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0}, "pT bin limits for pions"};

  // Configurables for delta y selection
  struct : ConfigurableGroup {
    Configurable<int> nBinsY{"nBinsY", 20, "Number of bins in y axis"};
    Configurable<int> nBinsDeltaY{"nBinsDeltaY", 20, "Number of bins in deltay axis"};
    Configurable<float> cfgYAcceptance{"cfgYAcceptance", 0.5f, "Rapidity acceptance"};
    Configurable<std::vector<float>> cfgDeltaYAcceptanceBins{"cfgDeltaYAcceptanceBins", {0.5f}, "Rapidity acceptance bins"};
  } yConfigs;

  // Configurable to apply efficiency online
  Configurable<bool> applyEfficiency{"applyEfficiency", false, "Use efficiency for filling histograms"};

  // Configurable for event mixing
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};

  // Configurables for CCDB
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository to use"};
  Configurable<std::string> ccdbEfficiencyPath{"ccdbEfficiencyPath", "Users/s/scannito/Efficiencies", "Correction path to file"};

  // Constants
  double massKa = o2::constants::physics::MassKPlus;
  double massPi = o2::constants::physics::MassPiPlus;
  double massK0S = o2::constants::physics::MassK0Short;
  double massLambda = o2::constants::physics::MassLambda0;

  // Filter on phi selected collisions
  Filter collisionFilter = aod::lf_selection_phi_collision::phimesonSel == true;

  // Defining filters on V0s (cannot filter on dynamic columns)
  Filter v0PreFilter = (nabs(aod::v0data::dcapostopv) > v0Configs.v0SettingDCAPosToPV && nabs(aod::v0data::dcanegtopv) > v0Configs.v0SettingDCANegToPV && aod::v0data::dcaV0daughters < v0Configs.v0SettingDCAV0Dau);

  // Defining the type of the collisions for data and MC
  using SelCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults, aod::PhimesonSelection>>;
  using SimCollisions = soa::Join<SelCollisions, aod::McCollisionLabels>;

  // Defining the type of the V0s and corresponding daughter tracks for data and MC
  using FullV0s = soa::Filtered<aod::V0Datas>;
  using FullMCV0s = soa::Join<FullV0s, aod::McV0Labels>;

  using V0DauTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullPi>;
  using V0DauMCTracks = soa::Join<V0DauTracks, aod::McTrackLabels>;

  // Defining the type of the tracks for data and MC
  using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTOFFullPi, aod::pidTOFFullKa>;
  using FilteredTracks = soa::Filtered<FullTracks>;

  using FullMCTracks = soa::Join<FullTracks, aod::McTrackLabels>;
  using FilteredMCTracks = soa::Filtered<FullMCTracks>;

  Partition<FullTracks> posTracks = aod::track::signed1Pt > trackConfigs.cfgCutCharge;
  Partition<FullTracks> negTracks = aod::track::signed1Pt < trackConfigs.cfgCutCharge;

  Partition<FullMCTracks> posMCTracks = aod::track::signed1Pt > trackConfigs.cfgCutCharge;
  Partition<FullMCTracks> negMCTracks = aod::track::signed1Pt < trackConfigs.cfgCutCharge;

  // Manual slicing
  Preslice<aod::Tracks> trackPerCollision = aod::track::collisionId;
  SliceCache cache;

  // Efficiency maps
  std::shared_ptr<TH3> effMapPhi = nullptr;
  std::shared_ptr<TH3> effMapK0S = nullptr;
  std::shared_ptr<TH3> effMapPionTPC = nullptr;
  std::shared_ptr<TH3> effMapPionTPCTOF = nullptr;

  void init(InitContext&)
  {
    AxisSpec deltayAxis = {yConfigs.nBinsDeltaY, -1.0f, 1.0f, "#Delta#it{y}"};
    AxisSpec deltaphiAxis = {72, -o2::constants::math::PIHalf, o2::constants::math::PIHalf * 3, "#Delta#varphi"};
    AxisSpec multAxis = {120, 0.0f, 120.0f, "centFT0M"};
    AxisSpec binnedmultAxis{(std::vector<double>)binsMult, "centFT0M"};
    AxisSpec massPhiAxis = {200, 0.9f, 1.2f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec pTPhiAxis = {120, 0.0f, 12.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTPhiAxis{(std::vector<double>)binspTPhi, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pTK0SAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTK0SAxis{(std::vector<double>)binspTK0S, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pTPiAxis = {50, 0.0f, 5.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTPiAxis{(std::vector<double>)binspTPi, "#it{p}_{T} (GeV/#it{c})"};

    histos.add("phi/h3PhiData", "Invariant mass of Phi in Data", kTH3F, {binnedmultAxis, binnedpTPhiAxis, massPhiAxis});
    histos.add("phiK0S/h5PhiK0SData2PartCorr", "Deltay vs deltaphi for Phi and K0Short in Data", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTK0SAxis, deltayAxis, deltaphiAxis});
    histos.add("phiPi/h5PhiPiData2PartCorr", "Deltay vs deltaphi for Phi and Pion in Data", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTPiAxis, deltayAxis, deltaphiAxis});

    // histos.add("phiK0S/h5PhiK0SDataNewProc", "2D Invariant mass of Phi and K0Short in Data", kTHnSparseF, {deltayAxis, binnedmultAxis, binnedpTK0SAxis, massK0SAxis, massPhiAxis});
    // histos.add("phiPi/h5PhiPiTPCDataNewProc", "Phi Invariant mass vs Pion nSigma TPC in Data", kTHnSparseF, {deltayAxis, binnedmultAxis, binnedpTPiAxis, nSigmaPiAxis, massPhiAxis});
    // histos.add("phiPi/h5PhiPiTOFDataNewProc", "Phi Invariant mass vs Pion nSigma TOF in Data", kTHnSparseF, {deltayAxis, binnedmultAxis, binnedpTPiAxis, nSigmaPiAxis, massPhiAxis});
  }

  // Compute weight based on efficiencies
  template <typename... BoundEffMaps>
  float computeWeight(const BoundEffMaps&... boundEffMaps)
  {
    if (!applyEfficiency)
      return 1.0f;

    float totalEfficiency = (boundEffMaps.getEfficiency() * ...);

    return totalEfficiency <= 0.0f ? 1.0f : 1.0f / totalEfficiency;
  }

  // Topological track selection for Phi's daughter tracks
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

  // Single track selection for strangeness sector
  template <typename T>
  bool selectionTrackStrangeness(const T& track)
  {
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsFound() < trackConfigs.minTPCnClsFound)
      return false;
    if (track.tpcNClsCrossedRows() < trackConfigs.minNCrossedRowsTPC)
      return false;
    if (track.tpcChi2NCl() > trackConfigs.maxChi2TPC)
      return false;

    if (std::abs(track.eta()) > trackConfigs.etaMax)
      return false;
    return true;
  }

  // V0 selection
  template <bool isMC, typename T1, typename T2>
  bool selectionV0(const T1& v0, const T2& collision)
  {
    const auto& posDaughterTrack = v0.template posTrack_as<V0DauTracks>();
    const auto& negDaughterTrack = v0.template negTrack_as<V0DauTracks>();

    if (!selectionTrackStrangeness(posDaughterTrack) || !selectionTrackStrangeness(negDaughterTrack))
      return false;

    if constexpr (!isMC) {
      if (std::abs(posDaughterTrack.tpcNSigmaPi()) > trackConfigs.nSigmaCutTPCPion)
        return false;
      if (std::abs(negDaughterTrack.tpcNSigmaPi()) > trackConfigs.nSigmaCutTPCPion)
        return false;
    }

    if (v0.v0cosPA() < v0Configs.v0SettingCosPA)
      return false;
    if (v0.v0radius() < v0Configs.v0SettingRadius)
      return false;
    if (v0.pt() < v0Configs.v0SettingMinPt)
      return false;

    if (v0Configs.cfgFurtherV0Selection) {
      if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massK0S > v0Configs.ctauK0s)
        return false;
      if (v0.qtarm() < (v0Configs.paramArmenterosCut * std::abs(v0.alpha())))
        return false;
      if (std::abs(v0.mLambda() - massLambda) < v0Configs.v0rejK0s)
        return false;
    }

    return true;
  }

  // Topological selection for pions
  template <typename T>
  bool selectionPion(const T& track)
  {
    if (!track.isGlobalTrackWoDCA())
      return false;

    if (track.itsNCls() < trackConfigs.minITSnCls)
      return false;
    if (track.tpcNClsFound() < trackConfigs.minTPCnClsFound)
      return false;

    if (track.pt() < trackConfigs.cMinPionPtcut)
      return false;

    if (trackConfigs.cfgIsTOFChecked && track.pt() >= trackConfigs.pTToUseTOF && !track.hasTOF())
      return false;

    if (std::abs(track.dcaXY()) > trackConfigs.cMaxDCArToPVPion->at(0) + (trackConfigs.cMaxDCArToPVPion->at(1) / std::pow(track.pt(), trackConfigs.cMaxDCArToPVPion->at(2))))
      return false;
    if (trackConfigs.cfgIsDCAzParameterized) {
      if (std::abs(track.dcaZ()) > trackConfigs.cMaxDCAzToPVPion->at(0) + (trackConfigs.cMaxDCAzToPVPion->at(1) / std::pow(track.pt(), trackConfigs.cMaxDCAzToPVPion->at(2))))
        return false;
    } else {
      if (std::abs(track.dcaZ()) > trackConfigs.cMaxDCAzToPVcut)
        return false;
    }
    return true;
  }

  void processPhiK0SPionDeltayDeltaphiData2D(SelCollisions::iterator const& collision, FullTracks const& fullTracks, FullV0s const& V0s, V0DauTracks const&)
  {
    float multiplicity = collision.centFT0M();

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    // Loop over all positive tracks
    for (const auto& track1 : posThisColl) {
      if (!selectionTrackResonance(track1) || !selectionPIDKaonpTdependent(track1))
        continue; // topological and PID selection

      // Loop over all negative tracks
      for (const auto& track2 : negThisColl) {
        if (!selectionTrackResonance(track2) || !selectionPIDKaonpTdependent(track2))
          continue; // topological and PID selection

        ROOT::Math::PxPyPzMVector recPhi = recMother(track1, track2, massKa, massKa);
        if (recPhi.Pt() < phiConfigs.minPhiPt)
          continue;
        if (recPhi.M() < phiConfigs.lowMPhi || recPhi.M() > phiConfigs.upMPhi)
          continue;
        if (std::abs(recPhi.Rapidity()) > yConfigs.cfgYAcceptance)
          continue;

        float weightPhi = computeWeight(BoundEfficiencyMap(effMapPhi, multiplicity, recPhi.Pt(), recPhi.Rapidity()));

        histos.fill(HIST("phi/h3PhiData"), multiplicity, recPhi.Pt(), recPhi.M(), weightPhi);

        // V0 already reconstructed by the builder
        for (const auto& v0 : V0s) {
          // Cut on V0 dynamic columns
          if (!selectionV0<false>(v0, collision))
            continue;
          if (std::abs(v0.yK0Short()) > yConfigs.cfgYAcceptance)
            continue;

          float weightPhiK0S = computeWeight(BoundEfficiencyMap(effMapPhi, multiplicity, recPhi.Pt(), recPhi.Rapidity()), BoundEfficiencyMap(effMapK0S, multiplicity, v0.pt(), v0.yK0Short()));

          histos.fill(HIST("phiK0S/h5PhiK0SData2PartCorr"), multiplicity, recPhi.Pt(), v0.pt(), recPhi.Rapidity() - v0.yK0Short(), recPhi.Phi() - v0.phi(), weightPhiK0S);
        }

        // Loop over all primary pion candidates
        for (const auto& track : fullTracks) {
          if (!selectionPion(track))
            continue;
          if (std::abs(track.rapidity(massPi)) > yConfigs.cfgYAcceptance)
            continue;

          float weightPhiPion = computeWeight(BoundEfficiencyMap(effMapPhi, multiplicity, recPhi.Pt(), recPhi.Rapidity()), track.pt() < trackConfigs.pTToUseTOF ? BoundEfficiencyMap(effMapPionTPC, multiplicity, track.pt(), track.rapidity(massPi)) : BoundEfficiencyMap(effMapPionTPCTOF, multiplicity, track.pt(), track.rapidity(massPi)));

          histos.fill(HIST("phiPi/h5PhiPiData2PartCorr"), multiplicity, recPhi.Pt(), track.pt(), recPhi.Rapidity() - track.rapidity(massPi), recPhi.Phi() - track.phi(), weightPhiPion);
        }
      }
    }
  }

  PROCESS_SWITCH(PhiStrangenessCorrelation, processPhiK0SPionDeltayDeltaphiData2D, "Process function for Phi-K0S and Phi-Pion Deltay and Deltaphi 2D Correlations in Data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PhiStrangenessCorrelation>(cfgc)};
}
