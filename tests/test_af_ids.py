import pytest
from apfid import parse_apfid, Apfid


class TestAlphafoldApfids:
    def test_af(self):
        apfid = parse_apfid("AF-A1KHV0-F1-V4_A")
        assert apfid.source == 'Alphafold'
        assert apfid.af_id is not None
        assert apfid.af_id.uniprot_id == 'A1KHV0'
        assert apfid.experiment_id == 'AF-A1KHV0-F1-model_v4'
        assert apfid.upper() == "AF-A1KHV0-F1-V4_A"

    def test_af_shorter(self):
        apfid = parse_apfid("AF-A1KHV0_A")
        assert apfid.source == 'Alphafold'
        assert apfid.af_id is not None
        assert apfid.af_id.uniprot_id == 'A1KHV0'
        assert apfid.experiment_id == 'AF-A1KHV0-F1-model_v4'
        assert apfid.upper() == "AF-A1KHV0-F1-V4_A"

    def test_af_older(self):
        apfid = parse_apfid("AF-A1KHV0-F1-V3_A")
        assert apfid.source == 'Alphafold'
        assert apfid.af_id is not None
        assert apfid.af_id.uniprot_id == 'A1KHV0'
        assert apfid.af_id.version == 3
        assert apfid.upper() == "AF-A1KHV0-F1-V3_A"

    def test_af_v2(self):
        apfid = parse_apfid("AF-A1KHV0_A20-50")
        assert apfid.source == 'Alphafold'
        assert apfid.start == 20
        assert apfid.end == 50
        assert apfid.af_id.uniprot_id == 'A1KHV0'
        assert apfid.version == 2