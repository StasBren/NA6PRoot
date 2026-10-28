#!/usr/bin/env python3
"""No Garfield dependency. Fake installation directories test env resolution."""
from pathlib import Path
import tempfile

from pi_garfield_runtime_env import garfield_runtime_env


def fake_database(prefix):
    folder=prefix/'share'/'Heed'/'database'
    folder.mkdir(parents=True,exist_ok=True)
    (folder/'cross_section_fake.txt').write_text('data')
    return folder


def main():
    with tempfile.TemporaryDirectory() as td:
        root=Path(td)
        default_install=root/'na6p'/'install'/'garfieldpp'
        explicit_install=root/'other'
        empty_install=root/'empty'
        default_db=fake_database(default_install)
        explicit_db=fake_database(explicit_install)
        (empty_install/'share'/'Heed'/'database').mkdir(parents=True)

        result=garfield_runtime_env(environ={},home=root)
        assert result['HEED_DATABASE']==str(default_db.resolve())
        assert result['GARFIELD_INSTALL']==str(default_install.resolve())

        result=garfield_runtime_env(environ={'GARFIELD_INSTALL':str(explicit_install)},home=root)
        assert result['HEED_DATABASE']==str(explicit_db.resolve())

        result=garfield_runtime_env(environ={'HEED_DATABASE':str(explicit_db)},home=root)
        assert result['HEED_DATABASE']==str(explicit_db.resolve())

        result=garfield_runtime_env(environ={'HEED_DATABASE':str(root/'nonexistent')},home=root)
        assert result['HEED_DATABASE']==str(default_db.resolve())
        result=garfield_runtime_env(environ={'HEED_DATABASE':str(empty_install/'share'/'Heed'/'database')},home=root)
        assert result['HEED_DATABASE']==str(default_db.resolve())

        no_defaults=root/'no_defaults'
        no_defaults.mkdir()
        try:
            garfield_runtime_env(environ={},home=no_defaults)
        except RuntimeError as exc:
            assert 'TrackHeed database missing' in str(exc)
        else:
            raise AssertionError('Missing database should fail closed')
    print('PASS Garfield/Heed runtime env: direct Python fallback, priority, invalid/empty and fail-closed')


if __name__=='__main__':
    main()
