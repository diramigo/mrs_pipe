import os
import re
from nipype.interfaces.base import BaseInterface, BaseInterfaceInputSpec, TraitedSpec, traits, File
from fsl_mrs.core.nifti_mrs import NIFTI_MRS
from fsl_mrs.core.basis import Basis

class Base_fsl_mrs_Interface(BaseInterface):
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

def mrs_io_decorator(out_file):
    # wrap an fsl_mrs command between read and write functions
    def decorator(func):
        def wrapper(*args, **kwargs):
            from fsl_mrs.utils import mrs_io
            from fsl_mrs.core.nifti_mrs import NIFTI_MRS

            if isinstance(args, tuple):
                args = list(args)
            else:
                args = [args]
            for index, arg in enumerate(args):
                args[index] = mrs_io.read_FID(arg)

            nifti_mrs = func(*args, **kwargs)
            
            if isinstance(nifti_mrs, tuple) and isinstance(out_file, list):
                for n_mrs, o_file in zip(nifti_mrs, out_file):
                    n_mrs.save(o_file)
            elif isinstance(nifti_mrs, NIFTI_MRS) and isinstance(out_file, str):
                nifti_mrs.save(out_file)
            else:
                raise Exception()

            return out_file
        
        return wrapper
    
    return decorator

        

# General Input and Output Specs

class Base_NIFTI_MRS_InputSpec(BaseInterfaceInputSpec):
    in_file = File(
        exists=True,
        desc="NIFTI_MRS data.", 
        mandatory=True
        )
    out_file = File(
        desc="NIFTI_MRS data.",
        mandatory=False
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

# ------------------ CoilCombine ------------------
class CoilCombineInputSpec(BaseInterfaceInputSpec):
    in_file = File(
        exists=True,
        desc="NIFTI_MRS data.", 
        mandatory=True
        )
    out_file = File(
        desc="NIFTI_MRS data.",
        mandatory=False
        )
    reference = File(
        desc="NIFTI_MRS water reference data.",
        mandatory=False,
        default=None
    )
    noise = traits.List(
        traits.Float,
        default=None,
        desc="Supply noise (NCoils x M) to estimate coil covariance (overridden by no_prewhiten)",
    )
    covariance = traits.List(
        traits.Float,
        default=None,
        desc="Supply coil-covariance for prewhitening (overridden by noise or no_prewhiten)",
    )
    no_prewiten = traits.Bool(
        desc="True to disable prewhitening."
        defualt=False,
        mandatory=False
    )
    report = traits.Either(
        traits.Bool, traits.Directory,
        desc="Provide output location as path to generate report. If set to True, uses working directory.",
        default=False,
        mandatory=False
    )
    report_all = traits.Bool(
        desc="True to output all indicies.",
        defualt=False,
        mandatory=False
        )

class CoilCombineOutputSpec(TraitedSpec):
    out_file = File(
        exists=True,
        desc="NIFTI_MRS data.",
        mandatory=True
        )
    report = traits.Directory(
        exists=True,
        desc="Preprocessing report.",
        mandatory=False
    )
    
class CoilCombine(Base_fsl_mrs_Interface):
    # Input and output specs
    input_spec = CoilCombineInputSpec
    output_spec = CoilCombineOutputSpec
    
    INTERFACE_NAME='coilcombine'

    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = super()._generate_out_file_name(self.inputs.in_file, self.inputs.out_file, self.INTERFACE_NAME)
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        # report
        if self.inputs.report:
            self.inputs.report = os.path.abspath(self.inputs.report)

        from fsl_mrs.utils.preproc import nifti_mrs_proc

        # run function
        self.inputs.out_file = mrs_io_decorator(self.inputs.out_file)(nifti_mrs_proc.coilcombine)(
            # mandatory file_names
            self.inputs.in_file, 
            # optional file names
            reference=self.inputs.reference,
            # optional parameters
            noise=self.inputs.noise,
            covariance=self.inputs.covariance,
            no_prewhiten=self.inputs.no_prewhiten,
            report=self.inputs.report,
            report_all=self.inputs.report_all
            )
        
        return runtime

    def _list_outputs(self):
        # get outputs
        outputs = self._outputs().get()
        outputs['out_file'] = self.inputs.out_file
        outputs['report'] = self.inputs.report
        return outputs


class AverageInputSpec(BaseInterfaceInputSpec):
    in_file = File(
        exists=True,
        desc="NIFTI_MRS data.", 
        mandatory=True
        )
    out_file = File(
        desc="NIFTI_MRS data.",
        mandatory=False
        )
    dim = traits.Enum(
        "DIM_DYN", "DIM_EDIT", 'all',
        desc="NIFTI-MRS dimension tag, or 'all'.",
        mandatory=True
    )
    report = traits.Either(
        traits.Bool, traits.Directory,
        desc="Provide output location as path to generate report. If set to True, uses working directory.",
        default=False,
        mandatory=False
    )
    report_all = traits.Bool(
        desc="True to output all indicies.",
        defualt=False,
        mandatory=False
        )

class AverageOutputSpec(TraitedSpec):
    out_file = File(
        exists=True,
        desc="NIFTI_MRS data.",
        mandatory=True
        )
    report = traits.Directory(
        exists=True,
        desc="Preprocessing report.",
        mandatory=False
    )

class Average(Base_fsl_mrs_Interface):
    # Input and output specs
    input_spec = AverageInputSpec
    output_spec = AverageOutputSpec
    
    INTERFACE_NAME='average'

    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = super()._generate_out_file_name(self.inputs.in_file, self.inputs.out_file, self.INTERFACE_NAME)
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        # report
        if self.inputs.report:
            self.inputs.report = os.path.abspath(self.inputs.report)

        from fsl_mrs.utils.preproc import nifti_mrs_proc

        # run function
        self.inputs.out_file = mrs_io_decorator(self.inputs.out_file)(nifti_mrs_proc.average)(
            # mandatory file_names
            self.inputs.in_file, 
            # mandatory parameters
            dim=self.inputs.dim
            # optional parameters
            # TODO: Add optional parameters
            )
        
        return runtime

    def _list_outputs(self):
        # get outputs
        outputs = self._outputs().get()
        outputs['out_file'] = self.inputs.out_file
        outputs['report'] = self.inputs.report
        return outputs


class AlignInputSpec(Base_fsl_mrs_Interface):
    in_file = File(
        exists=True,
        desc="NIFTI_MRS data.", 
        mandatory=True
        )
    out_file = File(
        desc="NIFTI_MRS data.",
        mandatory=False
        )
    dim = traits.Enum(
        "DIM_DYN", "DIM_EDIT", 'all',
        desc="NIFTI-MRS dimension tag, or 'all'.",
        mandatory=True
        )
    window = traits.Int(
        desc="Window size.",
        mandatory=False
        default=None
        )
    # target: I am not sure how this argument is supposed to work.
    ppmlim = traits.Tuple(
        (traits.Float, traits.Float), 
        desc="ppm search limits.", 
        mandatory=False,
        default=None
        )
    niter = traits.Int(
        desc="niter: Number of total iterations", 
        default=2, 
        mandatory=False
        )
    report = traits.Either(
        traits.Bool, traits.Directory,
        desc="Provide output location as path to generate report. If set to True, uses working directory.",
        default=False,
        mandatory=False
    )
    report_all = traits.Bool(
        desc="True to output all indicies.",
        defualt=False,
        mandatory=False
        )

class AlignOutputSpec(TraitedSpec):
    out_file = File(
        exists=True,
        desc="NIFTI_MRS data.",
        mandatory=True
        )
    report = traits.Directory(
        exists=True,
        desc="Preprocessing report.",
        mandatory=False
    )


class Align(Base_fsl_mrs_Interface):
    # Input and output specs
    input_spec = AlignInputSpec
    output_spec = AlignOutputSpec
    
    INTERFACE_NAME='align'

    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = super()._generate_out_file_name(self.inputs.in_file, self.inputs.out_file, self.INTERFACE_NAME)
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        # report
        if self.inputs.report:
            self.inputs.report = os.path.abspath(self.inputs.report)

        from fsl_mrs.utils.preproc import nifti_mrs_proc

        # run function
        self.inputs.out_file = mrs_io_decorator(self.inputs.out_file)(nifti_mrs_proc.align)(
            # mandatory file_names
            self.inputs.in_file, 
            # mandatory parameters
            self.inputs.dim,
            # optional parameters
            window=self.inputs.window,
            ppmlim=self.inputs.ppmlim,
            niter=self.inputs.niter,
            report=self.inputs.report,
            report_all=self.inputs.report_all
            )
        
        return runtime

    def _list_outputs(self):
        # get outputs
        outputs = self._outputs().get()
        outputs['out_file'] = self.inputs.out_file
        outputs['report'] = self.inputs.report
        return outputs


class EddyCurrentCorrectionInputSpec(Base_NIFTI_MRS_InputSpec, Ref_NIFTI_MRS_InputSpec):
    pass


class EddyCurrentCorrection(Base_fsl_mrs_Interface):
    input_spec = EddyCurrentCorrectionInputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    INTERFACE_NAME='ecc'

    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = super()._generate_out_file_name(self.inputs.in_file, self.inputs.out_file, self.INTERFACE_NAME)
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        from fsl_mrs.utils.preproc import nifti_mrs_proc

        # run function
        self.inputs.out_file = mrs_io_decorator(self.inputs.out_file)(nifti_mrs_proc.ecc)(
            # mandatory file_names
            self.inputs.in_file, 
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


class PhaseCorrect(Base_fsl_mrs_Interface):
    input_spec = PhaseCorrectInputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    INTERFACE_NAME='phasecorrrect'

    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = super()._generate_out_file_name(self.inputs.in_file, self.inputs.out_file, self.INTERFACE_NAME)
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        from fsl_mrs.utils.preproc import nifti_mrs_proc

        # run function
        self.inputs.out_file = mrs_io_decorator(self.inputs.out_file)(nifti_mrs_proc.phase_correct)(
            # mandatory file_names
            self.inputs.in_file, 
            # mandatory parameters
            ppmlim=self.inputs.ppmlim
            # optional parameters
            # TODO: Add optional parameters
            )
        
        return runtime

    def _list_outputs(self):
        # get outputs
        outputs = self._outputs().get()
        outputs['out_file'] = self.inputs.out_file
        return outputs
    

class PhaseCorrect_Creatine_ppmlim(Base_fsl_mrs_Interface):
    input_spec = Base_NIFTI_MRS_InputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    INTERFACE_NAME='phasecorrectCr'

    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = super()._generate_out_file_name(self.inputs.in_file, self.inputs.out_file, self.INTERFACE_NAME)
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        from fsl_mrs.utils.preproc import nifti_mrs_proc

        # run function
        self.inputs.out_file = mrs_io_decorator(self.inputs.out_file)(nifti_mrs_proc.phase_correct)(
            # mandatory file_names
            self.inputs.in_file, 
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
    ppm_ref = traits.Float(
        desc="Reference shift that peak will be moved to.", 
        mandatory=True
        )
    peak_search = traits.Tuple(
        (traits.Float,traits.Float), 
        desc="Search for peak between these ppm limits e.g. (2.8, 3.2) for tCr.", 
        mandatory=True
        )

class ShiftToReferenceInputSpec(Base_NIFTI_MRS_InputSpec, ppm_ref_InputSpec):
    pass


class ShiftToReference(Base_fsl_mrs_Interface):
    input_spec = ShiftToReferenceInputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    INTERFACE_NAME='shift2reference'

    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = super()._generate_out_file_name(self.inputs.in_file, self.inputs.out_file, self.INTERFACE_NAME)
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        from fsl_mrs.utils.preproc import nifti_mrs_proc

        # run function
        self.inputs.out_file = mrs_io_decorator(self.inputs.out_file)(nifti_mrs_proc.shift_to_reference)(
            # mandatory file_names
            self.inputs.in_file, 
            # mandatory parameters
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
    

class ShiftToCreatine(Base_fsl_mrs_Interface):
    input_spec = Base_NIFTI_MRS_InputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    INTERFACE_NAME='shift2Cr'

    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = super()._generate_out_file_name(self.inputs.in_file, self.inputs.out_file, self.INTERFACE_NAME)
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        from fsl_mrs.utils.preproc import nifti_mrs_proc

        # run function
        self.inputs.out_file = mrs_io_decorator(self.inputs.out_file)(nifti_mrs_proc.shift_to_reference)(
            # mandatory file_names
            self.inputs.in_file, 
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
    limits = traits.Tuple(
        (traits.Float(), traits.Float()), 
         desc="ppm limits between which peaks will be removed." , 
         mandatory=True
         )
    limit_units = traits.Enum(
        'ppm', 'ppm+shift', 'Hz', 
        desc='units of ppm_limits.', 
        mandatory=False
        )

class base_removewater_InputSpec(BaseInterfaceInputSpec):
    limits = traits.Tuple(
        -.25, .25, 
        desc="ppm limits between which peaks will be removed." , 
        usedefault=True, 
        mandatory=False
        )
    limit_units = traits.Enum(
        'ppm', 'ppm+shift', 'Hz', 
        desc='units of ppm_limits.', 
        default='ppm', 
        usedefault=True, 
        mandatory=False
        )

class RemovePeakInputSpec(Base_NIFTI_MRS_InputSpec, base_removepeak_InputSpec):
    pass

class RemoveWaterInputSpec(Base_NIFTI_MRS_InputSpec, base_removewater_InputSpec):
    pass


class RemovePeaks(Base_fsl_mrs_Interface):
    input_spec = RemovePeakInputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    INTERFACE_NAME='removepeaks'

    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = super()._generate_out_file_name(self.inputs.in_file, self.inputs.out_file, self.INTERFACE_NAME)
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        from fsl_mrs.utils.preproc import nifti_mrs_proc

        # run function
        self.inputs.out_file = mrs_io_decorator(self.inputs.out_file)(nifti_mrs_proc.remove_peaks)(
            # mandatory file_names
            self.inputs.in_file, 
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

class RemoveWater(Base_fsl_mrs_Interface):
    input_spec = RemoveWaterInputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    INTERFACE_NAME='removeH2O'

    def __init__(self, **inputs):
        super().__init__(**inputs)
        os.environ["OMP_NUM_THREADS"] = "4" # export OMP_NUM_THREADS=4
        os.environ["OPENBLAS_NUM_THREADS"] = "4" # export OPENBLAS_NUM_THREADS=4 
        os.environ["MKL_NUM_THREADS"] = "6" # export MKL_NUM_THREADS=6
        os.environ["VECLIB_MAXIMUM_THREADS"] = "4" # export VECLIB_MAXIMUM_THREADS=4
        os.environ["NUMEXPR_NUM_THREADS"] = "6" # export NUMEXPR_NUM_THREADS=6

    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = super()._generate_out_file_name(self.inputs.in_file, self.inputs.out_file, self.INTERFACE_NAME)
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        from fsl_mrs.utils.preproc import nifti_mrs_proc


        # run function
        self.inputs.out_file = mrs_io_decorator(self.inputs.out_file)(nifti_mrs_proc.remove_peaks)(
            # mandatory file_names
            self.inputs.in_file, 
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
    

class AlignByDynamicFit(Base_fsl_mrs_Interface):

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

def hermes_ref_remove_zero_mean_transients(ref):
    """
    Function for removing zero_mean transients in TwinsMX HERMES data.
    """
    
    from fsl_mrs.utils import mrs_io
    from fsl_mrs.utils.preproc import nifti_mrs_proc
    import fsl_mrs.core.nifti_mrs as nifti_mrs_tools
    import numpy as np


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

    return ref


class _HERMESRefRemoveZeroMeanTransients(Base_fsl_mrs_Interface):
    """
    For TwinsMX HERMES ref data. Can't guarantee that this is a necesarry or useful step for other HERMES data.
    """
    input_spec = Base_NIFTI_MRS_InputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    INTERFACE_NAME='rm0meanTransients'

    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = super()._generate_out_file_name(self.inputs.in_file, self.inputs.out_file, self.INTERFACE_NAME)
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        from fsl_mrs.utils.preproc import nifti_mrs_proc

        # run function
        self.inputs.out_file = mrs_io_decorator(self.inputs.out_file)(hermes_ref_remove_zero_mean_transients)(
            # mandatory file_names
            self.inputs.in_file, 
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

class SplitEditSubspectra(Base_fsl_mrs_Interface):
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



class HERMESAlignYEdit(Base_fsl_mrs_Interface):
    # Input and output specs
    input_spec = Base_NIFTI_MRS_InputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    INTERFACE_NAME='hermesYalign'

    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = super()._generate_out_file_name(self.inputs.in_file, self.inputs.out_file, self.INTERFACE_NAME)
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        from fsl_mrs.utils.preproc import nifti_mrs_proc

        # run function
        self.inputs.out_file = mrs_io_decorator(self.inputs.out_file)(hermes_align_y_edit)(
            # mandatory file_names
            self.inputs.in_file
            )
        
        return runtime

    def _list_outputs(self):
        # get outputs
        outputs = self._outputs().get()
        outputs['out_file'] = self.inputs.out_file
        return outputs


def hermes_sort_subspectra(nifti_mrs):
    """
    Sort HERMES DIM_EDIT spectra.

    ecc_ref is optional, but must be preprocessed.
    """
    import pandas as pd
    from fsl_mrs.utils.preproc import nifti_mrs_proc
    from fsl_mrs.core import nifti_mrs as nifti_mrs_tools
    from fsl_mrs.utils import mrs_io
    import numpy as np

    nifti_mrs_preproc = nifti_mrs.copy()

    if 'DIM_DYN' in nifti_mrs.dim_tags:
        averaged = False
    else:
        averaged = True
    
    # align DIM_DYN
    if not averaged:
        nifti_mrs_preproc = nifti_mrs_proc.align(nifti_mrs_preproc, dim='DIM_DYN')
        nifti_mrs_preproc = nifti_mrs_proc.average(nifti_mrs_preproc, 'DIM_DYN')
    
    # phase correction
    nifti_mrs_preproc = nifti_mrs_proc.phase_correct(nifti_mrs_preproc, (2.9,3.1))

    # yshift
    nifti_mrs_preproc = hermes_align_y_edit(nifti_mrs_preproc)

    # Identify subspectra type (GABA_on = low NAA signal, GSH_on = low Asp signal)
    peak_naa = []
    peak_h20 = []
    for subspectrum in split_edit_subspectra(nifti_mrs_preproc):
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

    # separate spectra
    subspectra_dict = {}
    for index, subspectra in enumerate(split_edit_subspectra(nifti_mrs)):
        subspectra_dict[peak_df.loc[index, 'pulse']] = subspectra

    # merge spectra
    nifti_mrs = nifti_mrs_tools.merge(
        (subspectra_dict['both_off'], subspectra_dict['gaba_on'], subspectra_dict['gsh_on'], subspectra_dict['both_on']), 
        dimension='DIM_EDIT'
        )
    
    return nifti_mrs


class hermes_sort_subspectra_InputSpec(Base_NIFTI_MRS_InputSpec, optional_ecc_NIFTI_MRS_InputSpec):
    pass

class HERMESSortSubspectra(Base_fsl_mrs_Interface):
    # Input and output specs
    input_spec = hermes_sort_subspectra_InputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    INTERFACE_NAME='hermesSort'

    def _run_interface(self, runtime):

        # output to tmp directory
        self.inputs.out_file = super()._generate_out_file_name(self.inputs.in_file, self.inputs.out_file, self.INTERFACE_NAME)
        self.inputs.out_file = os.path.abspath(self.inputs.out_file)

        from fsl_mrs.utils.preproc import nifti_mrs_proc

        # run function
        self.inputs.out_file = mrs_io_decorator(self.inputs.out_file)(hermes_sort_subspectra)(
            # mandatory file_names
            self.inputs.in_file, 
            )
        
        return runtime

    def _list_outputs(self):
        # get outputs
        outputs = self._outputs().get()
        outputs['out_file'] = self.inputs.out_file
        return outputs


class hermes_sum_InputSpec(BaseInterfaceInputSpec):
    in_file = File(
        exists=True,
        desc="NIFTI_MRS data.", 
        mandatory=True
        )
    out_dir = traits.Directory(
        exists=True, 
        desc="Output directory.",
        default = None
        )

class hermes_sum_OutputSpec(TraitedSpec):
    summation = File(
        exists=True,
        desc="NIFTI_MRS data. Sum of all hermes subspectra."
    )
    gaba = File(
        exists=True,
        desc="NIFTI_MRS data. GABA On - GABA Off."
    )
    gsh = File(
        exists=True,
        desc="NIFTI_MRS data. GSH On - GSH Off."
    )

def hermes_edit_sum(nifti_mrs):
    """
    Get HERMES edit data.
    """
    from fsl_mrs.utils import mrs_io
    from fsl_mrs.utils.preproc import nifti_mrs_proc
    
    a,b,c,d = split_edit_subspectra(nifti_mrs)

    summation = nifti_mrs_proc.add(nifti_mrs_proc.add(a,b), nifti_mrs_proc.add(c,d))
    gaba = nifti_mrs_proc.subtract(nifti_mrs_proc.add(b, d), nifti_mrs_proc.add(c,a))
    gsh = nifti_mrs_proc.subtract(nifti_mrs_proc.add(c, d), nifti_mrs_proc.add(b,a))

    return summation, gaba, gsh


class HERMES_Edit_Sum(Base_fsl_mrs_Interface):
    # Input and output specs
    input_spec = hermes_sum_InputSpec
    output_spec = hermes_sum_OutputSpec
    
    INTERFACE_NAME='hermesSum'

    def _run_interface(self, runtime):

        if not isinstance(self.inputs.out_dir, str):
            self.inputs.out_dir = os.getcwd()

        # output to tmp directory
        self._summation = super()._generate_out_file_name(self.inputs.in_file, None, "sum")
        self._summation  = os.path.abspath(os.path.join(self.inputs.out_dir, self._summation))

        self._gaba = super()._generate_out_file_name(self.inputs.in_file, None, "gaba")
        self._gaba  = os.path.abspath(os.path.join(self.inputs.out_dir, self._gaba))

        self._gsh = super()._generate_out_file_name(self.inputs.in_file, None, "gsh")
        self._gsh  = os.path.abspath(os.path.join(self.inputs.out_dir, self._gsh))

        output_list = [self._summation, self._gaba, self._gsh]

        # run function
        self._summation, self._gaba, self._gsh = mrs_io_decorator(output_list)(hermes_edit_sum)(
            # mandatory file_names
            self.inputs.in_file, 
            )
        
        return runtime

    def _list_outputs(self):
        # get outputs
        outputs = self._outputs().get()
        outputs['summation'] = self._summation
        outputs['gaba'] = self._gaba
        outputs['gsh'] = self._gsh
        return outputs
    

