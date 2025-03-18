import pytest
from apfid import parse_apfid, Apfid

class TestApfidParsing:
    """Tests for the parse_apfid function"""
    
    def test_long_v1(self):
        apfid = parse_apfid('1YSI_A111_A191')
        assert apfid.version == 1
        assert apfid.experiment_id == '1YSI'
        assert apfid.chain_id == 'A'
        assert apfid.start == 111
        assert apfid.end == 191
        assert apfid.lower() == '1ysi_A111_A191'

    def test_short_v1(self):
        apfid = parse_apfid('1YSI_A')
        assert apfid.version == 1
        assert apfid.experiment_id == '1YSI'
        assert apfid.chain_id == 'A'
        assert apfid.start is None
        assert apfid.end is None
        assert apfid.lower() == '1ysi_A'

    def test_v2(self):
        apfid_lines = [
            '1YSI:10_A111_191',
            '1YSI:10_A111-191',
            '1YSI:10_A',
        ]
        for line in apfid_lines:
            apfid = parse_apfid(line)
            assert apfid.version == 2
            assert apfid.experiment_id == '1YSI'
            assert apfid.chain_id == 'A'
            assert apfid.model == 10

    def test_source(self):
        assert parse_apfid('1YSI_A111_A191').source == 'PDB'
        assert parse_apfid('USR-123QWE_A').source == 'UserUpload'
        assert parse_apfid('AF-123QWE_A').source == 'Alphafold'

    def test_scratch(self):
        apfid = Apfid(
            experiment_id='1YSI',
            chain_id='A',
            start=111,
            end=191,
            version=1
        )
        assert apfid.upper() == '1YSI_A111_A191'
        apfid.set_version(2)
        assert apfid.version == 2
        assert apfid.upper() == '1YSI_A111_191'
