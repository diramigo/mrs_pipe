import os
import re
from nipype.interfaces.base import BaseInterface


class BaseMRSInterface(BaseInterface):
    # BIDS pattern to name files
    BIDS_KEY_VALUES = r"(sub-[a-zA-Z0-9]+)(_ses-[a-zA-Z0-9]+)?(_task-[a-zA-Z0-9]+)?(_acq-[a-zA-Z0-9]+)?(_nuc-[a-zA-Z0-9]+)?(_voi-[a-zA-Z0-9]+)?(_rec-[a-zA-Z0-9]+)?(_run-[a-zA-Z0-9]+)?(_echo-[a-zA-Z0-9]+)?(_inv-[a-zA-Z0-9]+)?"
    BIDS_DERIVATIVE =  r"(_desc-[a-zA-Z0-9]+)?"
    BIDS_SUFFIX_EXTENSION = r"(_[a-zA-Z0-9]+\.nii)(\.gz)?"
    BIDS_PATTERN = BIDS_KEY_VALUES+BIDS_DERIVATIVE+BIDS_SUFFIX_EXTENSION

    def _is_bids_valid(self, filename):
        # identify if filename is bids valid
        return bool(re.match(self.BIDS_PATTERN, os.path.basename(filename)))
    
    def _generate_bids_derivative_name(self, filename, desc):
        # generate a bids derivative name
        bids_key_values = "".join(re.findall(self.BIDS_KEY_VALUES, filename)[-1])
        bids_suffix_extension = "".join(re.findall(self.BIDS_SUFFIX_EXTENSION, filename)[-1])
        bids_desc = f'_desc-{desc}'
        return bids_key_values + bids_desc + bids_suffix_extension
    
    def _generate_out_file_name(self, in_file, out_file, interface_name):
        # given parameters, choose a name for output 
        out_bids_compliant = False
        out_is_str = False
        if isinstance(out_file, str):
            out_is_str = True
            if re.match('^[a-zA-Z0-9]+$', out_file):
                out_bids_compliant = True

        if self._is_bids_valid(in_file) and out_bids_compliant:
            return self._generate_bids_derivative_name(in_file, out_file)
        
        elif self._is_bids_valid(in_file) and not out_is_str and interface_name:
           return self._generate_bids_derivative_name(in_file, interface_name)
        
        elif not self._is_bids_valid(in_file) and not out_is_str and interface_name:
            return interface_name + '.nii.gz'
        else:
            return out_file
        
