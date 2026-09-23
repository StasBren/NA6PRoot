#!/bin/bash

# The variables below contain the number NEVENTS_BCK of background events in the sample and the number of signal events NEVENTS_SIG that one wants to merge with each background event. N.B.: the background simulation must be run in advance and its outputs (kinematics, VT hits, MS hits) renamed to XXX_bck.root

NEVENTS_BCK=10
NEVENTS_SIG=30

for ((i=0; i<NEVENTS_BCK; i++)); do
  echo "Running na6psim with background event $i"

  export BKG_EVENT=$i

  na6psim -n $NEVENTS_SIG \
          -g genDimuonBgEvent.C \
          -V embed.C \
	  -r -1 \
	  -v 0

  root -l -b -q "MergeHits_embed.C+($i)"	  
	  
done

OUTFILESH=""
OUTFILESK=""
for ((i=0; i<NEVENTS_BCK; i++)); do
  OUTFILESK+=" MCKine_mix_${i}.root"
  OUTFILESHM+=" HitsMuonSpec_mix_${i}.root"
  OUTFILESHV+=" HitsVerTel_mix_${i}.root"
done

# ---- MCKine file ----
if [ -f MCKine_mix_all.root ]; then
  read -p "MCKine_mix_all.root exists. Delete and recreate? [y/N] " ans
  if [[ $ans =~ ^[Yy]$ ]]; then
    rm MCKine_mix_all.root
    hadd -j 4 MCKine_mix_all.root $OUTFILESK
  else
    echo "Skipping MCKine merge."
  fi
else
  hadd -j 4 MCKine_mix_all.root $OUTFILESK
fi

# ---- Hits file muons----
if [ -f HitsMuonSpec_mix_all.root ]; then
  read -p "HitsMuonSpec_mix_all.root exists. Delete and recreate? [y/N] " ans
  if [[ $ans =~ ^[Yy]$ ]]; then
    rm HitsMuonSpec_mix_all.root
    hadd -j 4 HitsMuonSpec_mix_all.root $OUTFILESHM
  else
    echo "Skipping Muon Hits merge."
  fi
else
  hadd -j 4 HitsMuonSpec_mix_all.root $OUTFILESHM
fi

# ---- Hits file vertex telescope ----
if [ -f HitsVerTel_mix_all.root ]; then
  read -p "HitsVerTel_mix_all.root exists. Delete and recreate? [y/N] " ans
  if [[ $ans =~ ^[Yy]$ ]]; then
    rm HitsVerTel_mix_all.root
    hadd -j 4 HitsVerTel_mix_all.root $OUTFILESHV
  else
    echo "Skipping Vertex Telescope Hits merge."
  fi
else
  hadd -j 4 HitsVerTel_mix_all.root $OUTFILESHV
fi
