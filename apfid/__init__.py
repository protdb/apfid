import logging
import re
from typing import Union
import warnings

_RG_APFID_V1_FULL = re.compile(
    r"(?P<experiment_id>[A-Za-z0-9-]+)_(?P<chain_id>[A-Za-z]{1,2})(?P<start>\d+)[-_](?P<chain2_id>[A-Zaz]{1,2})(?P<end>\d+)"
)
_RG_APFID_V1_CHAIN = re.compile(
    r"(?P<experiment_id>[A-Za-z0-9-]+)_(?P<chain_id>[A-Za-z]{1,2})"
)
_RG_APFID_V2_CHAIN = re.compile(
    r"(?P<experiment_id>[A-Za-z0-9-]+):?(?P<model>\d*)_(?P<chain_id>[A-Za-z]{1,2})"
)
_RG_APFID_V2_FULL = re.compile(
    r"(?P<experiment_id>[A-Za-z0-9-]+):?(?P<model>\d*)_(?P<chain_id>[A-Za-z]{1,2})_?(?P<start>\d+)[-_](?P<chain2_id>[A-Za-z]{0,2})(?P<end>\d+)"
)

class Apfid:
    apfid: str
    version: int
    experiment_id: str
    chain_id: str
    chain2_id: str | None
    model: int = 0
    start: Union[int, None]
    end: Union[int, None]
    id_type: str
    source: str

    def __init__(self, experiment_id='', chain_id='', start=None, end=None, model=0, chain2_id=None, apfid=None, version=1):
        logging.debug(f"{experiment_id} {chain_id} {start} {end}")
        if apfid is None:
            self.experiment_id = experiment_id
            self.chain_id = chain_id
            self.chain2_id = chain2_id
            self.model = model
            if model is not None and model > 0:
                self.version = 2
            else:
                self.version = version
            if start == end:
                self.start = None
                self.end = None
            else:
                self.start = start
                self.end = end
            self.apfid = self._make_apfid()
        else:
            warnings.warn(DeprecationWarning('passing apfid is deprecated, use parse_apfid() instead'))
            self.apfid = apfid
            self._parse_apfid()
        self._set_experiment_type()

    def set_version(self, version):
        self.version = version
        self._make_apfid()

    def _set_experiment_type(self):
        if len(self.experiment_id) == 4:
            self.source = "PDB"
        else:
            if self.experiment_id.startswith('AF'):
                self.source = "Alphafold"
            elif self.experiment_id.startswith('USR'):
                self.source = "UserUpload"
            else:
                self.source = "Unknown"
        self.id_type = self.source

    def _make_apfid(self, lower=False, version=None):
        if version is None:
            version = self.version
        exp_id = self.experiment_id.lower() if lower else self.experiment_id.upper()

        if version == 1:
            if self.start == self.end:
                return f"{exp_id}_{self.chain_id}"
            else:
                return f"{exp_id}_{self.chain_id}{self.start}_{self.chain_id}{self.end}"
        elif version == 2:
            apfid_str = f'{exp_id}'
            if self.model > 0:
                apfid_str += f':{self.model}'
            apfid_str += f'_{self.chain_id}'
            if self.start != self.end and self.start is not None and self.end is not None:
                apfid_str += f'{self.start}'
                if self.chain2_id != self.chain_id and self.chain2_id is not None:
                    apfid_str += f'_{self.chain2_id}{self.end}'
                else:
                    apfid_str += f'_{self.end}'
            return apfid_str
        else:
            raise ValueError(f"Unsupported version: {version}")

    def _parse_apfid(self):
        split = self.apfid.split("_")
        self.experiment_id = split[0]
        if len(split) == 2:
            self.chain_id = split[1]
            self.start = None
            self.end = None
        else:
            self.chain_id = split[1][0]
            self.start = int(split[1][1:])
            self.end = int(split[2][1:])

    def __str__(self):
        return self.apfid

    def upper(self):
        return self._make_apfid(lower=False)

    def lower(self):
        return self._make_apfid(lower=True)


def parse_apfid(apfid: str) -> Apfid:
    res = {
        'experiment_id': '',
        'model': 0,
        'start': None,
        'end': None,
        'chain_id': None,
        'chain2_id': None
    }
    for rg, v in ((_RG_APFID_V1_FULL, 1), (_RG_APFID_V2_FULL, 2), (_RG_APFID_V1_CHAIN, 1), (_RG_APFID_V2_CHAIN, 2)):
        match = rg.match(apfid)
        if match:
            res.update(match.groupdict())
            res['version'] = v
            break
    if not res['experiment_id'] or not res['chain_id']:
        raise ValueError(f"Invalid apfid {apfid}")
    if not res['model']:
        res['model'] = 0
    else:
        res['model'] = int(res['model'])
    res['start'] = int(res['start']) if res['start'] is not None else None
    res['end'] = int(res['end']) if res['end'] is not None else None
    return Apfid(**res)


if __name__ == '__main__':
    print(parse_apfid('1YSI_A111-191').upper())