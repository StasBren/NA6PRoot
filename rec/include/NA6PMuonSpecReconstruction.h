// NA6PCCopyright

// Based on:
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

#ifndef NA6P_MUONSPEC_RECONSTRUCTION_H
#define NA6P_MUONSPEC_RECONSTRUCTION_H

#include <string>
#include <Rtypes.h>

// Class to steer the VT reconstruction

#include "NA6PMuonSpecCluster.h"
#include "NA6PMCTruthContainer.h"
#include "NA6PMCComposedLabel.h"
#include "NA6PTrack.h"
#include "NA6PVertex.h"
#include "NA6PReconstruction.h"

class TFile;
class TTree;
class NA6PMuonSpecHit;
class NA6PMuonSpecVertexerTracklets;
class NA6PTrackerCA;

struct MatchCandidate {
  int track;
  int segment;
  float score;
};

class NA6PMuonSpecReconstruction : public NA6PReconstruction
{
 public:
  NA6PMuonSpecReconstruction();
  ~NA6PMuonSpecReconstruction() override;

  bool initTracker();
  // methods to steer cluster reconstruction
  void createClustersOutput() override;
  void clearClusters() override
  {
    mClusters.clear();
    mCluMCLabels.clear_andfreememory();
  }
  void writeClusters() override;
  void closeClustersOutput() override;
  std::vector<NA6PMuonSpecCluster>& getClusters() { return *hClusPtr; }
  const std::vector<NA6PMuonSpecCluster>& getClusters() const { return *hClusPtr; }
  void useOwnedClusterStorage() { hClusPtr = &mClusters; }

  // fast method to smear the hits bypassing digitization and cluster finder
  void setClusterSpaceResolutionX(double clures) { mCluResX = clures; }
  void setClusterSpaceResolutionY(double clures) { mCluResY = clures; }
  void hitsToRecPoints(const std::vector<NA6PMuonSpecHit>& hits, int evID = 0);
  NA6PTrackerCA* getTracker() const { return mMSTracker.get(); }

  // methods to steer tracking
  void setClusters(std::vector<NA6PMuonSpecCluster>& clusters);
  void setClustersMCLabels(NA6PMCTruthContainer& mcCluLabels)
  {
    hCluMCLabelsPtr = &mcCluLabels;
  }

  void createTracksOutput() override;
  void clearTracks() override
  {
    mTracks.clear();
    mTrkMCLabels.clear();
  }
  void writeTracks() override;
  void closeTracksOutput() override;
  void runTracking();
  void runMSTrackMIDTrackletMatching();

 private:
  std::vector<NA6PMuonSpecCluster> mClusters, *hClusPtr = &mClusters;              // vector of clusters
  NA6PMCTruthContainer mCluMCLabels, *hCluMCLabelsPtr = &mCluMCLabels;             // cluster MC labels
  TFile* mClusFile = nullptr;                                                      // file with clusters
  TTree* mClusTree = nullptr;                                                      // tree of clusters
  double mCluResX = 100.e-4;                                                       // cluster resolution, cm (for fast simu)
  double mCluResY = 500.e-4;                                                       // cluster resolution, cm (for fast simu)
  std::vector<NA6PTrack> mTracks, *hTrackPtr = &mTracks;                           // vector of tracks
  std::vector<NA6PMCComposedLabel> mTrkMCLabels, *hTrkMCLabelsPtr = &mTrkMCLabels; // track MC labels
  TFile* mTrackFile = nullptr;                                                     // file with tracks
  TTree* mTrackTree = nullptr;                                                     // tree of tracks
  std::unique_ptr<NA6PTrackerCA> mMSTracker;                                       // tracker

  ClassDefNV(NA6PMuonSpecReconstruction, 1);
};

#endif
