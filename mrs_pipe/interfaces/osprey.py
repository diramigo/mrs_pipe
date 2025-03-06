import os
import subprocess
from nipype.interfaces.base import BaseInterfaceInputSpec, TraitedSpec, traits, File
from mrs_pipe.interfaces.base import BaseMRSInterface

def _run_matlab_command(command):
    cmd = f"matlab -nodesktop -nosplash -batch \"run('{command}')\""
    print(cmd)
    os.system(cmd)
    return cmd

def _get_spectral_registration_script(
        in_file, 
        out_file, 
        osprey_path, 
        spm_path,
        script_path,
        method=["Probabilistic","Robust"], 
        sequence=["HERMES", "MEGA"]
        ):
    
    if method == 'Probabilistic':
        method = 'op_probabSpecReg'
    elif method == 'Robust':
        method = 'op_robSpecReg'
    else:
        raise Exception("Not a valid spectral registration method.")

    
    matlab_command =  f"""
addpath(genpath('{osprey_path}'));
addpath(genpath('{spm_path}'));
raw = io_loadspec_niimrs('{in_file}');
seq='{sequence}';
ver_osp = ['Osprey ' getCurrentVersion().Version];
raw = osp_add_nii_mrs_field(raw,ver_osp);
[refShift_ind_ini]=op_preref(raw,seq);
[raw, fs, phs, weights, driftPre, driftPost] = {method}(raw, seq, 0,refShift_ind_ini);
io_writeniimrs(raw, '{out_file}', DIM_EDIT=['DIM_DYN']);
    """

    with open(script_path, 'w') as script:
        script.write(matlab_command)
    
    return script_path
    
class SpectralRegistrationInputSpec(BaseInterfaceInputSpec):
    in_file = File(
        exists=True,
        desc="NIFTI_MRS data.", 
        mandatory=True
        )
    method = traits.Either(
        "Probabilistic", "Robust",
        desc="Method for spectral registration.",
        mandatory=True
    )
    sequence = traits.Either(
        "HERMES", "MEGA",
        desc="Sequence name.",
        mandatory=True
    )
    out_file = File(
        desc="NIFTI_MRS data.",
        mandatory=False
        )
    osprey_path = traits.Directory(
        desc="Osprey directory.",
        exists=True,
        mandatory=False
        )
    spm_path = traits.Directory(
        desc="SPM directory.",
        exists=True,
        mandatory=False
        )
    
class SpectralRegistrationOutputSpec(TraitedSpec):
    out_file = File(
        exists=True,
        desc="NIFTI_MRS data.",
        mandatory=True
        )

class SpectralRegistration(BaseMRSInterface):
    # Input and output specs
    input_spec = SpectralRegistrationInputSpec
    output_spec = SpectralRegistrationOutputSpec
    
    INTERFACE_NAME='specreg'

    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = super()._generate_out_file_name(self.inputs.in_file, self.inputs.out_file, self.INTERFACE_NAME)
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        # osprey path
        if not self.inputs.osprey_path:
            self.inputs.osprey_path = os.getenv('osprey')
            if not os.path.exists(self.inputs.osprey_path):
                raise FileExistsError("Osprey path does not exist.")

        # run function
        spectral_registration_script = _get_spectral_registration_script(
            self.inputs.in_file,
            self.inputs.out_file,
            self.inputs.osprey_path,
            self.inputs.spm_path,
            os.path.join(os.getcwd(), 'script.m'),
            self.inputs.method,
            self.inputs.sequence
            )
        
        x = _run_matlab_command(spectral_registration_script)

        return runtime

    def _list_outputs(self):
        # get outputs
        outputs = self._outputs().get()
        outputs['out_file'] = self.inputs.out_file
        return outputs
