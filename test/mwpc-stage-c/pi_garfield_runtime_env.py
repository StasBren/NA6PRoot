#!/usr/bin/env python3
"""Resolve Garfield++/Heed runtime for subprocesses, also from direct Python.

This is a pure standard-library helper. It never modifies the login shell,
installs packages or changes the user's source/physics output directories.
"""
from __future__ import annotations

import os
from pathlib import Path


def garfield_runtime_env(*, environ=None, home=None):
    """Return a validated copy of the environment for TrackHeed subprocesses.

    Source-of-truth order: explicit valid HEED_DATABASE, valid
    GARFIELD_INSTALL/share/Heed/database, GARFIELD_HOME install/source tree,
    then the established ~/na6p/install/garfieldpp project install prefix.
    Never accept an existing directory that contains no database files.
    """
    env = dict(os.environ if environ is None else environ)
    user_home = Path(home).expanduser() if home is not None else Path.home()
    pref = user_home / 'na6p' / 'install' / 'garfieldpp'
    candidates = []

    def add(path, install=None):
        if not path:
            return
        p = Path(path).expanduser()
        if (p, install) not in candidates:
            candidates.append((p, install))

    add(env.get('HEED_DATABASE'), env.get('GARFIELD_INSTALL'))
    if env.get('GARFIELD_INSTALL'):
        install = Path(env['GARFIELD_INSTALL']).expanduser()
        add(install / 'share' / 'Heed' / 'database', install)
    if env.get('GARFIELD_HOME'):
        source = Path(env['GARFIELD_HOME']).expanduser()
        add(source / 'install' / 'share' / 'Heed' / 'database', source / 'install')
        add(source / 'share' / 'Heed' / 'database', source)
        add(source / 'Heed' / 'heed++' / 'database', source)
    add(pref / 'share' / 'Heed' / 'database', pref)
    # Also accept the local Garfield source database when a WSL install
    # has binaries but has not installed the Heed data into share/.
    source_default = user_home / 'na6p' / 'src' / 'garfieldpp'
    add(source_default / 'Heed' / 'heed++' / 'database', pref)

    for path, install in candidates:
        if not path.is_dir():
            continue
        try:
            has_regular_file = any(p.is_file() for p in path.rglob('*'))
        except OSError:
            continue
        if not has_regular_file:
            continue
        env['HEED_DATABASE'] = str(path.resolve())
        if install and Path(install).is_dir():
            env['GARFIELD_INSTALL'] = str(Path(install).resolve())
        return env

    inspected = '\n  '.join(str(p) for p, _ in candidates) or '(no paths)'
    raise RuntimeError(
        'TrackHeed database missing. Searched these locations:\n  '
        + inspected
        + '\nRun: find "$HOME/na6p" -type d -path "*/share/Heed/database" -print'
        + '\nIf the folder exists elsewhere, set HEED_DATABASE to that absolute directory.'
        + '\nDo not start B1 avalanches until this passes.'
    )


def heed_environment_check(*, environ=None, home=None):
    e = garfield_runtime_env(environ=environ, home=home)
    return f"PASS Garfield runtime: HEED_DATABASE={e['HEED_DATABASE']}"
