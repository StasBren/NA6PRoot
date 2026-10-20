#!/usr/bin/env python3
"""Synthetic test for true Geant4 local u-w midpoint/angle correlations."""
from __future__ import annotations
import csv,json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from pi_native_local_uw_density import process_csv,merge,build,NU,NW,GEOM,FIELDS,selections,resolve_root_executable,activate_root_for_export

def row(det,st,theta,u,w,pdg=13,path=6):
    return dict(event=0,hit_index=0,track_id=0,det_id=det,station=st,pdg=pdg,
      edep_eV=1000.,path_length_mm=path,
      u_in_mm=u,u_out_mm=u,w_in_mm=w,w_out_mm=w,
      theta_u_in_deg=theta,theta_w_in_deg=0.,
      theta_u_out_deg=theta,theta_w_out_deg=0.)

def test():
    assert len(GEOM)==298 and NU==16 and NW==24
    assert GEOM[0].station==0 and GEOM[32].station==1
    assert dict(selections(0.))==dict(positive=True,far10=False,near2=True)
    assert dict(selections(10.))==dict(positive=True,far10=True,near2=False)
    assert not dict(selections(2.))["near2"]
    with TemporaryDirectory(prefix="pi_native_uw_") as tmp:
        root=Path(tmp);directory=root/"exports";directory.mkdir()
        rows=[row(0,0,0.,0.,0.),row(0,0,15.,-30.,50.),
          row(0,0,-20.,60.,-200.),row(0,0,1.,400.,0.),
          row(0,0,1.,0.,0.,path=0.),row(0,0,0.,0.,0.,pdg=211),
          row(32,1,-1.,0.,30.)]
        for channel in ("jpsi","phi","omega"):
            with (directory/f"{channel}_geant_gas_crossings.csv").open(
              "w",newline="",encoding="utf-8") as f:
                writer=csv.DictWriter(f,fieldnames=list(rows[0]))
                writer.writeheader();writer.writerows(rows)
        info,bad=process_csv(directory/"jpsi_geant_gas_crossings.csv")
        assert info["0"]["n"]==4 and info["0"]["outside"]==1
        assert info["32"]["n"]==1 and info["32"]["outside"]==0
        assert bad=={"invalid_native_hit":1,"outside_physical_module_or_nonfinite":1}
        assert sum(info["0"]["all"].values())==3
        assert sum(info["0"]["selected"]["theta_u_in_deg"]["positive"].values())==2
        assert sum(info["0"]["selected"]["theta_u_in_deg"]["far10"].values())==2
        assert sum(info["0"]["selected"]["theta_u_in_deg"]["near2"].values())==1
        combined={};merge(combined,info);merge(combined,info)
        assert combined["0"]["n"]==8
        assert sum(combined["0"]["all"].values())==6
        assert all(sum(combined["0"]["angles"][field])==8 for field in FIELDS)
        args=SimpleNamespace(mode="first",first=directory,out=root/"output",
            campaign=root)
        path=build(args)
        payload=json.loads(path.read_text(encoding="utf-8"))
        assert payload["meta"]["scope"]=="FIRST_CHUNK_3000_EVENTS_PER_SOURCE"
        assert all(len(m)==298 for m in payload["modules"].values())
        module=payload["modules"]["phi40-bore"]["0"]
        assert module["n"]==4 and module["outside"]==1
        assert sum(k for _,k in module["all"])==3
        assert sum(k for _,k in module["selected"]["theta_u_in_deg"]
                       ["near2"])==1
        # Regression for the WSL failure in full 201-chunk mode:
        # installed NA6PRoot is on PATH, but the root command is not.
        # Resolve real ROOT from the already-existing CMakeCache, without
        # assuming the interactive bash profile has sourced thisroot.sh.
        root_prefix=root/"existing-root-install"
        root_bin=root_prefix/"bin";root_bin.mkdir(parents=True)
        (root_prefix/"lib/cmake/ROOT").mkdir(parents=True)
        executable=root_bin/"root"
        executable.write_text("#!/bin/sh\nexit 0\n",encoding="utf-8")
        executable.chmod(0o755)
        cmake=root/"existing-cmake-build";cmake.mkdir()
        (cmake/"CMakeCache.txt").write_text(
            "ROOT_DIR:PATH="+str(root_prefix/"lib/cmake/ROOT")+"\n",
            encoding="utf-8")
        isolated={"PATH":"/does-not-contain-ROOT","LD_LIBRARY_PATH":"/other/lib"}
        resolved=resolve_root_executable(build_dirs=[cmake],env=isolated)
        assert resolved==executable.resolve()
        active=activate_root_for_export(build_dirs=[cmake],env=isolated)
        assert active==str(executable.resolve())
        assert isolated["ROOTSYS"]==str(root_prefix)
        assert isolated["PI_ROOT_EXECUTABLE"]==active
        assert isolated["PATH"].split(":")[0]==str(root_bin)
        assert str(root_prefix/"lib") in isolated["LD_LIBRARY_PATH"].split(":")
        assert "/other/lib" in isolated["LD_LIBRARY_PATH"].split(":")
        failed={**isolated,"PI_ROOT_EXECUTABLE":"/not/a/ROOT/bin/root"}
        try:
            resolve_root_executable(build_dirs=[cmake],env=failed)
        except FileNotFoundError as exc:
            assert "PI_ROOT_EXECUTABLE" in str(exc)
        else:
            raise AssertionError("Bad explicit ROOT path silently accepted")
        print("PASS: locate missing ROOT from existing CMake ROOT_DIR")
        print("PASS: activate rootcling PATH/ROOTSYS/shared libs in WSL")
        print("PASS: reject stale explicitly configured ROOT path")
        print("PASS: signed muon angle cuts, local gas aperture, native QA rejection")
        print("PASS: all-source 298-chamber 2D/1D midpoint exposure conservation")

if __name__=="__main__":test()
