import os

from nipype.interfaces.base import BaseInterface, BaseInterfaceInputSpec, TraitedSpec, traits, File
from fsl_mrs.core.nifti_mrs import NIFTI_MRS
from fsl_mrs.core.basis import Basis

# Functions

def read_NIFTI_MRS(path):
    """
    Read NIFTI MRS
    """
    
    from fsl_mrs.utils import mrs_io
    
    NIFTI_MRS = mrs_io.read_FID(path)
    
    return NIFTI_MRS


def read_basis(path, type=['press', 'hermes']):
    """
    Read NIFTI MRS
    """
    from pathlib import Path
    from fsl_mrs.utils import mrs_io
    
    match type:
        case 'press':
            basis = mrs_io.basis(path)
            return basis
        case 'hermes':
            a = mrs_io.basis(Path(path) / 'a')
            b = mrs_io.basis(Path(path) / 'b')
            c = mrs_io.basis(Path(path) / 'c')
            d = mrs_io.basis(Path(path) / 'd')
            return a,b,c,d
        case _:
            return None
        

# General Input and Output Specs

class Base_NIFTI_MRS_InputSpec(BaseInterfaceInputSpec):
    in_file = File(
        exists=True,
        desc="NIFTI_MRS data.", 
        mandatory=True
        )
    out_file = File(
        desc="NIFTI_MRS data.",
        mandatory=True
        )
    
class Base_NIFTI_MRS_OutputSpec(TraitedSpec):
    out_file = File(
        exists=True,
        desc="NIFTI_MRS data.",
        mandatory=True
        )
    
class Ref_NIFTI_MRS_InputSpec(BaseInterfaceInputSpec):
    ref = File(
        exists=True,
        desc="NIFTI_MRS water reference data.", 
        mandatory=True
        )
    
class optional_ecc_NIFTI_MRS_InputSpec(BaseInterfaceInputSpec):
    ecc_ref = File(
        exists=True,
        desc="NIFTI_MRS reference data for eddy current correction.", 
        mandatory=False
        )

class Dim_NIFTI_MRS_InputSpec(BaseInterfaceInputSpec):
    dim = traits.Enum(
        "DIM_DYN", "DIM_EDIT", 'all',
        desc="NIFTI-MRS dimension tag, or 'all'.",
        mandatory=True
    )
    
class ppm_NIFTI_MRS_InputSpec(BaseInterfaceInputSpec):
    ppmlim = traits.Tuple((traits.Float, traits.Float), desc="ppm range.", mandatory=True)

class optional_ppm_NIFTI_MRS_InputSpec(BaseInterfaceInputSpec):
    ppmlim = traits.Either(None, traits.Tuple((traits.Float, traits.Float)), desc="ppm range.", usedefault=True, mandatory=False)

class Base_Basis_MRS_InputSpec(BaseInterfaceInputSpec):
    basis = traits.Either(
        traits.Directory(exists=True), 
        traits.List(traits.Directory(exists=True)), 
        desc="Basis to fit, or list of basis matching the dimension to align.", 
        mandatory=True
        )
    
class Model_NIFTI_MRS_InputSpec(BaseInterfaceInputSpec):
    baseline_order = traits.Int(desc="Baseline order for modelling.", mandatory=True)
    model = traits.Enum(
        "voigt", "lorentzian",
        desc="Model type.",
        mandatory=True
    )



# Interfaces

# ------------------ Align ------------------
class AlignInputSpec(Base_NIFTI_MRS_InputSpec, Dim_NIFTI_MRS_InputSpec, optional_ppm_NIFTI_MRS_InputSpec):
    # TODO: Add optional parameters
    # Optional inputs
    window=traits.Int(desc="Window size.", default=None, mandatory=False)
    # target=None,
    # ppmlim=traits.Tuple((traits.Float(), traits.Float()), desc="ppm search limits.", default=None, mandatory=False)
    niter=traits.Int(desc="niter: Number of total iterations", default=2, mandatory=False)
    # figure=False,
    # report=None,
    # report_all=False,
    # :param NIFTI_MRS data: Data to align
    # :param str dim: NIFTI-MRS dimension tag, or 'all'
    # :param int window: Window size.
    # :param target: Optional target FID
    # :param ppmlim: ppm search limits.
    # :param int niter: Number of total iterations
    # :param figure: True to show figure.
    # :param report: Provide output location as path to generate report
    # :param report_all: True to output all indicies

def _align(in_file, out_file, dim, ppmlim):
    # TODO: Add optional parameters
    """
    Align NIFTI_MRS.
    """
    from fsl_mrs.utils.mrs_io import read_FID
    from fsl_mrs.utils.preproc import nifti_mrs_proc

    nifti_mrs = read_FID(in_file)
    
    if dim in nifti_mrs.dim_tags:
        print("Performing align in DIM_DYN.")
        nifti_mrs = nifti_mrs_proc.align(nifti_mrs, dim, ppmlim=ppmlim)

        print(f"Save NIfTI MRS at {out_file}.")
        nifti_mrs.save(out_file)
    else: 
        # if skipping step is ok
        print(f"No {dim}, skipping align step.")
        out_file = in_file

    return out_file

class Align(BaseInterface):
    # Input and output specs
    input_spec = AlignInputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        # run function
        self.inputs.out_file = _align(
            # mandatory file_names
            self.inputs.in_file, 
            self.inputs.out_file,
            # mandatory parameters
            self.inputs.dim,
            # optional parameters
            self.inputs.ppmlim
            # TODO: Add optional parameters
            )
        
        return runtime

    def _list_outputs(self):
        # get outputs
        outputs = self._outputs().get()
        outputs['out_file'] = self.inputs.out_file
        return outputs

# ------------------ Average ------------------

class AverageInputSpec(Base_NIFTI_MRS_InputSpec, Dim_NIFTI_MRS_InputSpec):
    pass

def _average(in_file, out_file, dim):
    """
    Average NIFTI_MRS.
    """
    from fsl_mrs.utils.preproc import nifti_mrs_proc
    from fsl_mrs.utils.mrs_io import read_FID

    nifti_mrs = read_FID(in_file)
    
    if dim in nifti_mrs.dim_tags:
        nifti_mrs = nifti_mrs_proc.average(nifti_mrs, dim)
        nifti_mrs.save(out_file)
    else:
        out_file = in_file

    return out_file

class Average(BaseInterface):
    # Input and output specs
    input_spec = AverageInputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        # run function
        self.inputs.out_file = _average(
            # mandatory file_names
            self.inputs.in_file, 
            self.inputs.out_file,
            # mandatory parameters
            self.inputs.dim
            # optional parameters
            # TODO: Add optional parameters
            )
        
        return runtime

    def _list_outputs(self):
        # get outputs
        outputs = self._outputs().get()
        outputs['out_file'] = self.inputs.out_file
        return outputs
    
# ------------------ EddyCurrentCorrection ------------------
class EddyCurrentCorrectionInputSpec(Base_NIFTI_MRS_InputSpec, Ref_NIFTI_MRS_InputSpec):
    pass

def _ecc(in_file, out_file, ecc_ref):
    """
    Perform eddy current correction.
    """
    from fsl_mrs.utils.preproc import nifti_mrs_proc
    from fsl_mrs.utils import mrs_io

    nifti_mrs = mrs_io.read_FID(in_file)
    ecc_ref = mrs_io.read_FID(ecc_ref)

    if ecc_ref is not None:
        nifti_mrs = nifti_mrs_proc.ecc(nifti_mrs, ecc_ref)
        nifti_mrs.save(out_file)
    else:
        out_file = in_file
    
    return out_file

class EddyCurrentCorrection(BaseInterface):
    input_spec = EddyCurrentCorrectionInputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        # run function
        self.inputs.out_file = _ecc(
            # mandatory file_names
            self.inputs.in_file, 
            self.inputs.out_file,
            # mandatory parameters
            self.inputs.ref
            # optional parameters
            # TODO: Add optional parameters
            )
        
        return runtime

    def _list_outputs(self):
        # get outputs
        outputs = self._outputs().get()
        outputs['out_file'] = self.inputs.out_file
        return outputs

# ------------------ PhaseCorrect ------------------
class PhaseCorrectInputSpec(Base_NIFTI_MRS_InputSpec, ppm_NIFTI_MRS_InputSpec):
    pass

def _phase_correct(in_file, out_file, ppmlim=(2.9,3.1)):
    """
    Phase correct NIFTI MRS.
    """
    from fsl_mrs.utils.preproc import nifti_mrs_proc
    from fsl_mrs.utils.mrs_io import read_FID

    nifti_mrs = read_FID(in_file)
    
    nifti_mrs = nifti_mrs_proc.phase_correct(nifti_mrs, ppmlim=ppmlim)

    nifti_mrs.save(out_file)
    
    return out_file

class PhaseCorrect(BaseInterface):
    input_spec = PhaseCorrectInputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        # run function
        self.inputs.out_file = _phase_correct(
            # mandatory file_names
            self.inputs.in_file, 
            self.inputs.out_file,
            # mandatory parameters
            self.inputs.ppmlim
            # optional parameters
            # TODO: Add optional parameters
            )
        
        return runtime

    def _list_outputs(self):
        # get outputs
        outputs = self._outputs().get()
        outputs['out_file'] = self.inputs.out_file
        return outputs
    

class PhaseCorrect_Creatine_ppmlim(BaseInterface):
    input_spec = Base_NIFTI_MRS_InputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        # run function
        self.inputs.out_file = _phase_correct(
            # mandatory file_names
            self.inputs.in_file, 
            self.inputs.out_file,
            # mandatory parameters
            ppmlim=(2.9,3.1)
            # optional parameters
            # TODO: Add optional parameters
            )
        
        return runtime

    def _list_outputs(self):
        # get outputs
        outputs = self._outputs().get()
        outputs['out_file'] = self.inputs.out_file
        return outputs
    
    
# ------------------ ShiftToReference ------------------
class ppm_ref_InputSpec(BaseInterfaceInputSpec):
    ppm_ref = traits.Float(desc="Reference shift that peak will be moved to.", mandatory=True)
    peak_search = traits.Tuple((traits.Float,traits.Float), desc="Search for peak between these ppm limits e.g. (2.8, 3.2) for tCr.", mandatory=True)

class ShiftToReferenceInputSpec(Base_NIFTI_MRS_InputSpec, ppm_ref_InputSpec):
    pass

def _shift_to_reference(in_file, out_file, ppm_ref, peak_search):
    """
    Shift to reference.
    """
    from fsl_mrs.utils.preproc import nifti_mrs_proc
    from fsl_mrs.utils import mrs_io

    nifti_mrs = mrs_io.read_FID(in_file)
    
    nifti_mrs = nifti_mrs_proc.shift_to_reference(
        nifti_mrs, 
        ppm_ref=ppm_ref, 
        peak_search=peak_search
        )
    
    nifti_mrs.save(out_file)

    return out_file

class ShiftToReference(BaseInterface):
    input_spec = ShiftToReferenceInputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        # run function
        self.inputs.out_file = _shift_to_reference(
            # mandatory file_names
            self.inputs.in_file, 
            self.inputs.out_file,
            # fixed parameters
            ppm_ref=self.inputs.ppm_ref,
            peak_search=self.inputs.peak_search,
            # optional parameters
            # TODO: Add optional parameters
            )
        
        return runtime

    def _list_outputs(self):
        # get outputs
        outputs = self._outputs().get()
        outputs['out_file'] = self.inputs.out_file
        return outputs
    

class ShiftToCreatine(BaseInterface):
    input_spec = Base_NIFTI_MRS_InputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        # run function
        self.inputs.out_file = _shift_to_reference(
            # mandatory file_names
            self.inputs.in_file, 
            self.inputs.out_file,
            # fixed parameters
            ppm_ref=3.027,
            peak_search=(2.9, 3.1),
            # optional parameters
            # TODO: Add optional parameters
            )
        
        return runtime

    def _list_outputs(self):
        # get outputs
        outputs = self._outputs().get()
        outputs['out_file'] = self.inputs.out_file
        return outputs

# ------------------ ShiftToReference ------------------
class base_removepeak_InputSpec(BaseInterfaceInputSpec):
    limits = traits.Tuple((traits.Float(), traits.Float()), desc="ppm limits between which peaks will be removed." , mandatory=True)
    limit_units = traits.Enum('ppm', 'ppm+shift', 'Hz', desc='units of ppm_limits.', mandatory=False)

class base_removewater_InputSpec(BaseInterfaceInputSpec):
    limits = traits.Tuple(-.25, .25, desc="ppm limits between which peaks will be removed." , usedefault=True, mandatory=False)
    limit_units = traits.Enum('ppm', 'ppm+shift', 'Hz', desc='units of ppm_limits.', default='ppm', usedefault=True, mandatory=False)

class RemovePeakInputSpec(Base_NIFTI_MRS_InputSpec, base_removepeak_InputSpec):
    pass

class RemoveWaterInputSpec(Base_NIFTI_MRS_InputSpec, base_removewater_InputSpec):
    pass

def _remove_peaks(in_file, out_file, limits, limit_units):
    """
    Remove peak.
    """
    from fsl_mrs.utils.preproc import nifti_mrs_proc
    from fsl_mrs.utils import mrs_io
    
    nifti_mrs = mrs_io.read_FID(in_file)

    nifti_mrs = nifti_mrs_proc.remove_peaks(nifti_mrs, limits=limits, limit_units=limit_units)
    
    nifti_mrs.save(out_file)

    return out_file

class RemovePeaks(BaseInterface):
    input_spec = RemovePeakInputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        # run function
        self.inputs.out_file = _remove_peaks(
            # mandatory file_names
            self.inputs.in_file, 
            self.inputs.out_file,
            # mandatory parameters
            limits=self.inputs.limits,
            limit_units=self.inputs.limit_units,
            # optional parameters
            # TODO: Add optional parameters
            )
        
        return runtime

    def _list_outputs(self):
        # get outputs
        outputs = self._outputs().get()
        outputs['out_file'] = self.inputs.out_file
        return outputs

class RemoveWater(BaseInterface):
    input_spec = RemoveWaterInputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        # run function
        self.inputs.out_file = _remove_peaks(
            # mandatory file_names
            self.inputs.in_file, 
            self.inputs.out_file,
            # fixed parameters
            limits=self.inputs.limits,
            limit_units=self.inputs.limit_units,
            # optional parameters
            # TODO: Add optional parameters
            )
        
        return runtime

    def _list_outputs(self):
        # get outputs
        outputs = self._outputs().get()
        outputs['out_file'] = self.inputs.out_file
        return outputs

# ------------------ DynamicAlign ------------------

class AlignByDynamicFit_InputSpec(Base_NIFTI_MRS_InputSpec, Base_Basis_MRS_InputSpec, Model_NIFTI_MRS_InputSpec, ppm_NIFTI_MRS_InputSpec):
    pass

class eps_phi_OutputSpec(TraitedSpec):
    eps = traits.Float(desc="Frequency shift")
    phi = traits.Float(desc="Zero-order phase.")

class AlignByDynamicFit_OutputSpec(Base_NIFTI_MRS_OutputSpec, eps_phi_OutputSpec):
    pass

def _dynamic_align(in_file:str, out_file: str, basis: list, baseline_order:int, model, ppmlim:tuple):
    """
    Perform eddy current correction.
    """
    from fsl_mrs.utils.preproc import dyn_based_proc as dproc
    from fsl_mrs.utils import mrs_io

    nifti_mrs = mrs_io.read_FID(in_file)

    # read basis
    if isinstance(basis, list):
        basis = [mrs_io.read_basis(single_basis) for single_basis in basis]
    else:
        basis = mrs_io.read_basis(basis)

    fitargs = {
        "baseline_order":baseline_order, 
        "model":model, 
        "ppmlim":ppmlim
        }

    nifti_mrs, eps, phi = dproc.align_by_dynamic_fit(nifti_mrs, basis, fitargs)
    
    nifti_mrs.save(out_file)

    return out_file, eps, phi
    

class AlignByDynamicFit(BaseInterface):

    input_spec = AlignByDynamicFit_InputSpec
    output_spec = AlignByDynamicFit_OutputSpec
    
    def _run_interface(self, runtime):
        # output to tmp directory
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        # run function
        self.inputs.out_file, self._eps, self._phi = _dynamic_align(
            # mandatory file_names
            self.inputs.in_file, 
            self.inputs.out_file,
            # mandatory directories / directory
            self.inputs.basis,
            # mandatory parameters
            baseline_order=self.inputs.baseline_order,
            model=self.inputs.model,
            ppmlim=self.inputs.ppmlim
            # optional parameters
            # TODO: Add optional parameters
            )
        
        return runtime

    def _list_outputs(self):
        # get outputs
        outputs = self._outputs().get()
        outputs['out_file'] = self.inputs.out_file
        outputs['eps'] = self._eps
        outputs['phi'] = self._phi
        return outputs


###################
###### Utils ######
###################

def _split_ref_dyn(ref):
    """
    Function to split ref DIM_DYN. Used to remove zero-mean transients in TwinsMX HERMES data.
    """
    import fsl_mrs.core.nifti_mrs as nifti_mrs_tools

    ref_list = []
    for _ref_index in range(ref.shape[4] - 1):
        ref_subspectra, ref = nifti_mrs_tools.split(ref, 'DIM_DYN', 0)
        ref_list.append(ref_subspectra)
    ref_list.append(ref)
    return tuple(ref_list)

def _hermes_ref_remove_zero_mean_transients(in_file, out_file):
    """
    Function for removing zero_mean transients in TwinsMX HERMES data.
    """
    
    from fsl_mrs.utils import mrs_io
    from fsl_mrs.utils.preproc import nifti_mrs_proc
    import fsl_mrs.core.nifti_mrs as nifti_mrs_tools
    import numpy as np

    # read 
    ref = mrs_io.read_FID(in_file)
    # only the first edit dimension has data on it
    ref, _empty = nifti_mrs_tools.split(ref, dimension='DIM_EDIT', index_or_indices=0)
    # remove dim_edit
    ref = nifti_mrs_proc.average(ref, 'DIM_EDIT')
    # get indices of non zero-mean transients
    indices = np.where(np.mean(ref[0,0,0,:, :].real, axis=0) != 0)[0]
    # extract non zero-mean transients
    ref_tuple = [subref for index, subref in enumerate(_split_ref_dyn(ref)) if index in indices]
    # merge them back
    ref = nifti_mrs_tools.merge(ref_tuple, dimension='DIM_DYN')

    # save nifti_mrs
    ref.save(out_file)

    return out_file


class _HERMESRefRemoveZeroMeanTransients(BaseInterface):
    """
    For TwinsMX HERMES ref data. Can't guarantee that this is a necesarry or useful step for other HERMES data.
    """
    input_spec = Base_NIFTI_MRS_InputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        # run function
        self.inputs.out_file = _hermes_ref_remove_zero_mean_transients(
            # mandatory file_names
            self.inputs.in_file, 
            self.inputs.out_file,
            )
        
        return runtime

    def _list_outputs(self):
        # get outputs
        outputs = self._outputs().get()
        outputs['out_file'] = self.inputs.out_file
        return outputs


def split_edit_subspectra(svs):
    """
    Split an SVS along the EDIT dimension.
    """
    from fsl_mrs.core import nifti_mrs as nifti_mrs_tools

    if "DIM_EDIT" not in svs.dim_tags:
        print("Not an edited sequence")
        return svs
    
    # first four dimensions are fixed, but DIM_EDIT could be the 5th,6th or 7th dimension
    for index, dim_tag in enumerate(svs.dim_tags):
        if dim_tag == "DIM_EDIT":
            break
    edit_dim_index = index + 4

    edit_list = []
    for _edit_index in range(svs.shape[edit_dim_index] - 1):
        sub_spectra, svs = nifti_mrs_tools.split(svs, 'DIM_EDIT', 0)
        edit_list.append(sub_spectra)

    edit_list.append(svs)

    return tuple(edit_list)

def _split_edit_subspectra(in_file, out_dir):
    """
    Split an SVS along the EDIT dimension. 
    """
    from fsl_mrs.utils import mrs_io

    extension = 'nii.gz'
    nifti_mrs = mrs_io.read_FID(in_file)

    subspectra = split_edit_subspectra(nifti_mrs)

    subspectra_filenames = []
    for subspectra_index in range(subspectra):
        subspectra_filenames.append(os.path.join(out_dir, subspectra_index + extension))
    
    for subspectrum, filename in zip(subspectra, subspectra_filenames):
        subspectrum.save(filename)

    return subspectra_filenames

class split_edit_subspectra_InputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, desc="NIFTI_MRS data.", mandatory=True)
    out_dir = traits.Directory(exists=True, desc="Output directory.", mandatory=False)
class split_edit_subspectra_OutputSpec(TraitedSpec):
    a = File(exists=True, desc="NIFTI_MRS data.", mandatory=True)
    b = File(exists=True, desc="NIFTI_MRS data.", mandatory=True)
    c = File(exists=False, desc="NIFTI_MRS data.", mandatory=False)
    d = File(exists=False, desc="NIFTI_MRS data.", mandatory=False)

class SplitEditSubspectra(BaseInterface):
    """
    Split DIM_EDIT.
    """
    input_spec = split_edit_subspectra_InputSpec
    output_spec = split_edit_subspectra_OutputSpec
    
    def _run_interface(self, runtime):

        # output to tmp directory
        if self.inputs.out_dir:
            self.inputs.out_dir = os.path.abspath(self.inputs.out_dir)
        else:
            self.inputs.out_dir = os.path.abspath(os.path.getcwd())


        # run function
        out_files = _split_edit_subspectra(
            # mandatory file_names
            self.inputs.in_file, 
            self.inputs.out_dir,
            )
        
        if len(out_files) >= 2:
            self._outputs.a = out_files[0]
            self._outputs.b = out_files[1]
        else:
            self._outputs.a = None
            self._outputs.b = None
        if len(out_files) == 4:
            self._outputs.c = out_files[2]
            self._outputs.d = out_files[3]
        else:
            self._outputs.c = None
            self._outputs.d = None
        
        return runtime

    def _list_outputs(self):
        # get outputs
        outputs = self._outputs().get()
        outputs['a'] = self.inputs.a
        outputs['b'] = self.inputs.b
        outputs['c'] = self.inputs.c
        outputs['d'] = self.inputs.d
        return outputs


def hermes_align_y_edit(svs):
    """
    Align subspectra in y-axis.
    """
    from fsl_mrs.utils.preproc import nifti_mrs_proc
    from fsl_mrs.utils.misc import SpecToFID
    from fsl_mrs.core import nifti_mrs as nifti_mrs_tools
    import numpy as np


    if 'DIM_DYN' in svs.dim_tags:
        averaged = False
    else:
        averaged = True
    
    # align DIM_DYN
    if not averaged:
        nifti_mrs = nifti_mrs_proc.average(svs, 'DIM_DYN')
    else:
        nifti_mrs = svs
    
    # Get mean spectrum of every subspectra
    spectrum_dict = {}
    for label, subspectra in zip(['A','B','C','D'], split_edit_subspectra(nifti_mrs)):
        spectrum_dict[label] = subspectra
        
    # get data in ppm-range of no interest
    spec_a = spectrum_dict['A'].mrs().get_spec([1, -20]).real
    spec_b = spectrum_dict['B'].mrs().get_spec([1, -20]).real
    spec_c = spectrum_dict['C'].mrs().get_spec([1, -20]).real
    spec_d = spectrum_dict['D'].mrs().get_spec([1, -20]).real

    # scaling factor
    median_a_b_diff = np.median(spec_a - spec_b)
    median_a_c_diff = np.median(spec_a - spec_c)
    median_a_d_diff = np.median(spec_a - spec_d)

    # scale non-prprocessed svs
    A, B, C, D = split_edit_subspectra(svs)
    
    if not averaged:
        for index, subspectrum in enumerate(B.mrs()):
            B[0,0,0,:,0,index] = SpecToFID(subspectrum.get_spec() + median_a_b_diff)
        
        for index, subspectrum in enumerate(C.mrs()):
            C[0,0,0,:,0,index] = SpecToFID(subspectrum.get_spec() + median_a_c_diff)
        
        for index, subspectrum in enumerate(D.mrs()):
            D[0,0,0,:,0,index] = SpecToFID(subspectrum.get_spec() + median_a_d_diff)
    else:
        B[0,0,0,:,0,0] = SpecToFID(B.mrs().get_spec() + median_a_b_diff)
        C[0,0,0,:,0,0] = SpecToFID(C.mrs().get_spec() + median_a_c_diff)
        D[0,0,0,:,0,0] = SpecToFID(D.mrs().get_spec() + median_a_d_diff)

    svs = nifti_mrs_tools.merge([A,B,C,D], 'DIM_EDIT')

    return svs

def _hermes_align_y_edit(in_file, out_file):
    """
    Align subspectra in y-axis.
    """
    from fsl_mrs.utils import mrs_io

    nifti_mrs = mrs_io.read_FID(in_file)

    nifti_mrs = hermes_align_y_edit(nifti_mrs)

    nifti_mrs.save(out_file)

    return out_file


class HERMESAlignYEdit(BaseInterface):
    # Input and output specs
    input_spec = Base_NIFTI_MRS_InputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        # run function
        self.inputs.out_file = _hermes_align_y_edit(
            # mandatory file_names
            self.inputs.in_file, 
            self.inputs.out_file,
            )
        
        return runtime

    def _list_outputs(self):
        # get outputs
        outputs = self._outputs().get()
        outputs['out_file'] = self.inputs.out_file
        return outputs


def _hermes_sort_subspectra(in_file, out_file, ecc_ref=None):
    """
    Sort HERMES DIM_EDIT spectra.

    ecc_ref is optional, but must be preprocessed.
    """
    import pandas as pd
    from fsl_mrs.utils.preproc import nifti_mrs_proc
    from fsl_mrs.core import nifti_mrs as nifti_mrs_tools
    from fsl_mrs.utils import mrs_io
    import numpy as np

    nifti_mrs = mrs_io.read_FID(in_file)

    if 'DIM_DYN' in nifti_mrs.dim_tags:
        averaged = False
    else:
        averaged = True
    
    # align DIM_DYN
    if not averaged:
        nifti_mrs = nifti_mrs_proc.align(nifti_mrs, dim='DIM_DYN')
        nifti_mrs = nifti_mrs_proc.average(nifti_mrs, 'DIM_DYN')
    
    # # ecc
    # if ecc_ref:
    #     ecc_ref = mrs_io.read_FID(ecc_ref)
    #     nifti_mrs = nifti_mrs_proc.ecc(nifti_mrs, ecc_ref)

    # # remove water
    # nifti_mrs = nifti_mrs_proc.remove_peaks(nifti_mrs, (-.25, .25), 'ppm')

    # # shift to creatine
    # nifti_mrs = nifti_mrs_proc.shift_to_reference(
    #     nifti_mrs, 
    #     ppm_ref=3.027, 
    #     peak_search=(2.9, 3.1)
    #     )
    
    # phase correction
    nifti_mrs = nifti_mrs_proc.phase_correct(nifti_mrs, (2.9,3.1))

    # yshift
    nifti_mrs = hermes_align_y_edit(nifti_mrs)

    # Identify subspectra type (GABA_on = low NAA signal, GSH_on = low Asp signal)
    peak_naa = []
    peak_h20 = []
    for subspectrum in split_edit_subspectra(nifti_mrs):
        # sum of NAA peak
        peak_naa.append(subspectrum.mrs().get_spec([2.2,1.8]).real.sum().real.max())
        # sum of Asp peak
        peak_h20.append(subspectrum.mrs().get_spec([5,4]).real.max())

    def which_pulse(x):
        if (not x['gsh_on']) and (not x['gaba_on']) :
            return 'both_off'
        elif (x['gsh_on']) and (x['gaba_on']):
            return 'both_on'
        elif (x['gsh_on']) and (not x['gaba_on']):
            return 'gsh_on'
        elif (not x['gsh_on']) and (x['gaba_on']):
            return 'gaba_on'
        else:
            return None

    peak_df = (pd.DataFrame({'subspectra':['a','b','c','d'], 'peak_naa':peak_naa, 'peak_h20':peak_h20})
                .assign(
                    gsh_on = lambda df_: df_['peak_h20'].rank() < 3,
                    gaba_on = lambda df_: df_.groupby('gsh_on')['peak_naa'].rank() < 2
                )
                .assign(pulse = lambda df_: df_.apply(which_pulse, axis=1))
                )

    #read file again
    nifti_mrs = mrs_io.read_FID(in_file)

    # separate spectra
    subspectra_dict = {}
    for index, subspectra in enumerate(split_edit_subspectra(nifti_mrs)):
        subspectra_dict[peak_df.loc[index, 'pulse']] = subspectra

    # merge spectra
    nifti_mrs = nifti_mrs_tools.merge(
        (subspectra_dict['both_off'], subspectra_dict['gaba_on'], subspectra_dict['gsh_on'], subspectra_dict['both_on']), 
        dimension='DIM_EDIT'
        )
    
    nifti_mrs.save(out_file)
    
    return out_file


class hermes_sort_subspectra_InputSpec(Base_NIFTI_MRS_InputSpec, optional_ecc_NIFTI_MRS_InputSpec):
    pass

class HERMESSortSubspectra(BaseInterface):
    # Input and output specs
    input_spec = hermes_sort_subspectra_InputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        # run function
        self.inputs.out_file = _hermes_sort_subspectra(
            # mandatory file_names
            self.inputs.in_file, 
            self.inputs.out_file,
            # optional file_names
            self.inputs.ecc_ref
            )
        
        return runtime

    def _list_outputs(self):
        # get outputs
        outputs = self._outputs().get()
        outputs['out_file'] = self.inputs.out_file
        return outputs
