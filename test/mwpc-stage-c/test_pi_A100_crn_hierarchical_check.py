#!/usr/bin/env python3
"""Synthetic offline checks for SHA-guarded A20/A100 CRN/bootstrap diagnostics.

Never launches Geant4, TrackHeed or Garfield physics and never writes original data.
"""
from pi_A100_crn_hierarchical_check import selftest

if __name__ == "__main__":
    selftest()
