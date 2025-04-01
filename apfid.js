/**
 * APFID JavaScript implementation
 * Ported from Python version
 */

// Regular expressions for parsing
const _RG_APFID_V1_FULL = /(?<experiment_id>[A-Za-z0-9-]+)_(?<chain_id>[A-Za-z]{1,2})(?<start>\d+)[-_](?<chain2_id>[A-Zaz]{1,2})(?<end>\d+)/;
const _RG_APFID_V1_CHAIN = /(?<experiment_id>[A-Za-z0-9-]+)_(?<chain_id>[A-Za-z]{1,2})/;
const _RG_APFID_V2_CHAIN = /(?<experiment_id>[A-Za-z0-9-]+):?(?<model>\d*)_(?<chain_id>[A-Za-z]{1,2})/;
const _RG_APFID_V2_FULL = /(?<experiment_id>[A-Za-z0-9-]+):?(?<model>\d*)_(?<chain_id>[A-Za-z]{1,2})_?(?<start>\d+)[-_](?<chain2_id>[A-Za-z]{0,2})(?<end>\d+)/;

/**
 * AlphafoldId class for handling Alphafold IDs
 */
class AlphafoldId {
  prefix = 'AF';
  uniprot_id = 'UNKNOWN';
  file_no = 1;
  version = 4;

  // Regular expressions for file number and version
  static _file_no_rg = /[Ff]{1}(?<no>\d+)/;
  static _version_rg = /[vV]{1}(?<no>\d+)/;

  /**
   * Create a new AlphafoldId instance
   * @param {string} alphafold_id - The Alphafold ID to parse
   */
  constructor(alphafold_id) {
    const af_params = alphafold_id.split(/[_-]/);
    if (af_params[0].toUpperCase() !== 'AF') {
      throw new Error(`Alphafold id ${alphafold_id} does not start with AF`);
    }
    
    this.uniprot_id = af_params[1];
    
    if (af_params.length > 2) {
      const fileMatch = AlphafoldId._file_no_rg.exec(af_params[2]);
      if (fileMatch && fileMatch.groups) {
        this.file_no = parseInt(fileMatch.groups.no);
      }
    }
    
    if (af_params.length > 3) {
      const versionMatch = AlphafoldId._version_rg.exec(af_params[af_params.length - 1]);
      if (versionMatch && versionMatch.groups) {
        this.version = parseInt(versionMatch.groups.no);
      }
    }
  }

  /**
   * Convert to string representation
   * @returns {string} String representation of the AlphafoldId
   */
  toString() {
    return `${this.prefix}-${this.uniprot_id}-F${this.file_no}-v${this.version}`;
  }

  /**
   * Get the download ID
   * @returns {string} The download ID
   */
  getDlId() {
    return `${this.prefix}-${this.uniprot_id}-F${this.file_no}-model_v${this.version}`;
  }

  /**
   * Get the PSSKB ID
   * @returns {string} The PSSKB ID
   */
  getPsskbId() {
    return `${this.prefix}-${this.uniprot_id}-F${this.file_no}-V${this.version}`;
  }
}

/**
 * Apfid class for handling protein fragment identifiers
 */
class Apfid {
  apfid = '';
  version = 1;
  experiment_id = '';
  chain_id = '';
  chain2_id = null;
  model = 0;
  start = null;
  end = null;
  id_type = '';
  source = '';
  af_id = null;

  /**
   * Create a new Apfid instance
   * @param {Object} options - Configuration options
   * @param {string} [options.experiment_id=''] - Experiment ID
   * @param {string} [options.chain_id=''] - Chain ID
   * @param {number|null} [options.start=null] - Start position
   * @param {number|null} [options.end=null] - End position
   * @param {number} [options.model=0] - Model number
   * @param {string|null} [options.chain2_id=null] - Second chain ID
   * @param {string|null} [options.apfid=null] - APFID string (deprecated, use parseApfid instead)
   * @param {number} [options.version=1] - APFID version
   */
  constructor({
    experiment_id = '',
    chain_id = '',
    start = null,
    end = null,
    model = 0,
    chain2_id = null,
    apfid = null,
    version = 1
  } = {}) {
    if (apfid === null) {
      this.experiment_id = experiment_id;
      this.chain_id = chain_id;
      this.chain2_id = chain2_id;
      this.model = model;
      
      if (model !== null && model > 0) {
        this.version = 2;
      } else {
        this.version = version;
      }
      
      if (start === end) {
        this.start = null;
        this.end = null;
      } else {
        this.start = start;
        this.end = end;
      }
      
      this.apfid = this._makeApfid();
    } else {
      console.warn('Warning: passing apfid is deprecated, use parseApfid() instead');
      this.apfid = apfid;
      this._parseApfid();
    }
    
    this._setExperimentType();
  }

  /**
   * Set the APFID version
   * @param {number} version - The version number
   */
  setVersion(version) {
    this.version = version;
    this._makeApfid();
  }

  /**
   * Determine the experiment type based on the experiment ID
   * @private
   */
  _setExperimentType() {
    if (this.experiment_id.length === 4) {
      this.source = "PDB";
    } else {
      if (this.experiment_id.startsWith('AF')) {
        this.source = "Alphafold";
        this.af_id = new AlphafoldId(this.experiment_id);
        this.experiment_id = this.af_id.getDlId();
      } else if (this.experiment_id.startsWith('USR')) {
        this.source = "UserUpload";
      } else {
        this.source = "Unknown";
      }
    }
    this.id_type = this.source; // compatibility
  }

  /**
   * Create an APFID string from the object properties
   * @param {boolean} [lower=false] - Convert to lowercase
   * @param {number|null} [version=null] - Override version
   * @returns {string} The APFID string
   * @private
   */
  _makeApfid(lower = false, version = null) {
    if (version === null) {
      version = this.version;
    }
    
    let exp_id;
    if (this.af_id !== null) {
      exp_id = this.af_id.getPsskbId();
    } else {
      exp_id = this.experiment_id;
    }
    
    exp_id = lower ? exp_id.toLowerCase() : exp_id.toUpperCase();

    if (version === 1) {
      if (this.start === this.end) {
        return `${exp_id}_${this.chain_id}`;
      } else {
        return `${exp_id}_${this.chain_id}${this.start}_${this.chain_id}${this.end}`;
      }
    } else if (version === 2) {
      let apfid_str = `${exp_id}`;
      
      if (this.model > 0) {
        apfid_str += `:${this.model}`;
      }
      
      apfid_str += `_${this.chain_id}`;
      
      if (this.start !== this.end && this.start !== null && this.end !== null) {
        apfid_str += `${this.start}`;
        if (this.chain2_id !== this.chain_id && this.chain2_id !== null) {
          apfid_str += `_${this.chain2_id}${this.end}`;
        } else {
          apfid_str += `_${this.end}`;
        }
      }
      
      return apfid_str;
    } else {
      throw new Error(`Unsupported version: ${version}`);
    }
  }

  /**
   * Parse an APFID string into object properties
   * @private
   */
  _parseApfid() {
    const split = this.apfid.split("_");
    this.experiment_id = split[0];
    
    if (split.length === 2) {
      this.chain_id = split[1];
      this.start = null;
      this.end = null;
    } else {
      this.chain_id = split[1][0];
      this.start = parseInt(split[1].substring(1));
      this.end = parseInt(split[2].substring(1));
    }
  }

  /**
   * Convert to string representation
   * @returns {string} The APFID string
   */
  toString() {
    return this.apfid;
  }

  /**
   * Get uppercase version of the APFID
   * @returns {string} Uppercase APFID
   */
  upper() {
    return this._makeApfid(false);
  }

  /**
   * Get lowercase version of the APFID
   * @returns {string} Lowercase APFID
   */
  lower() {
    return this._makeApfid(true);
  }
}

/**
 * Parse an APFID string and return an Apfid object
 * @param {string} apfid - The APFID string to parse
 * @returns {Apfid} The parsed Apfid object
 */
function parseApfid(apfid) {
  const res = {
    experiment_id: '',
    model: 0,
    start: null,
    end: null,
    chain_id: null,
    chain2_id: null
  };
  
  const regexVersionPairs = [
    [_RG_APFID_V1_FULL, 1],
    [_RG_APFID_V2_FULL, 2],
    [_RG_APFID_V1_CHAIN, 1],
    [_RG_APFID_V2_CHAIN, 2]
  ];

  let matched = false;
  
  for (const [regex, version] of regexVersionPairs) {
    const match = regex.exec(apfid);
    if (match && match.groups) {
      Object.assign(res, match.groups);
      res.version = version;
      matched = true;
      break;
    }
  }
  
  if (!res.experiment_id || !res.chain_id) {
    throw new Error(`Invalid apfid ${apfid}`);
  }
  
  if (!res.model) {
    res.model = 0;
  } else {
    res.model = parseInt(res.model);
  }
  
  res.start = res.start !== null ? parseInt(res.start) : null;
  res.end = res.end !== null ? parseInt(res.end) : null;
  
  return new Apfid(res);
}

// Export the classes and functions for use in other modules
export {
  AlphafoldId,
  Apfid,
  parseApfid
};

// Example usage (similar to the Python __main__ block)
if (typeof require !== 'undefined' && require.main === module) {
  console.log(parseApfid('1YSI_A111-191').upper());
}