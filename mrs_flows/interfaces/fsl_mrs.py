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
        

def preproc_spectrum(svs, ecc_ref=None):
    """
    Preprocess HERMES spectra that have been aligned within DIM_DYN.
    """
    if ecc_ref is not None:
        svs = NIFTI_MRS_proc.ecc(svs, ecc_ref)

    hlsvdlimits = [-0.25, 0.25]
    svs = NIFTI_MRS_proc.remove_peaks(svs, hlsvdlimits, limit_units='ppm')
    svs = NIFTI_MRS_proc.shift_to_reference(svs, ppm_ref=3.027, peak_search=(2.9, 3.1))
    svs = NIFTI_MRS_proc.phase_correct(svs, ppmlim=(2.9, 3.1))
    return svs


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

class Dim_NIFTI_MRS_InputSpec(BaseInterfaceInputSpec):
    dim = traits.Enum(
        "DIM_DYN", "DIM_EDIT", 'all',
        desc="NIFTI-MRS dimension tag, or 'all'.",
        mandatory=True
    )
    
class ppm_NIFTI_MRS_InputSpec(BaseInterfaceInputSpec):
    ppmlim = traits.Tuple((traits.Float, traits.Float), desc="ppm range.", mandatory=True)
    
class Model_NIFTI_MRS_InputSpec(BaseInterfaceInputSpec):
    basis = traits.Either(
        traits.Instance(Basis),
        traits.Tuple(traits.Instance(Basis)),
        desc="MRS Basis set.",
        mandatory=True
        )
    baseline_order = traits.Int(desc="Baseline order for modelling.", mandatory=True)
    model = traits.Enum(
        "voigt", "lorentzian",
        desc="Model type.",
        mandatory=True
    )

# Interfaces

# ------------------ Align ------------------
class AlignInputSpec(Base_NIFTI_MRS_InputSpec, Dim_NIFTI_MRS_InputSpec):
    # TODO: Add optional parameters
    # Optional inputs
    window=traits.Int(desc="Window size.", default=None, mandatory=False)
    # target=None,
    ppmlim=traits.Tuple((traits.Float(), traits.Float()), desc="ppm search limits.", default=None, mandatory=False)
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

def _align(in_file, out_file, dim='DIM_DYN'):
    # TODO: Add optional parameters
    """
    Align NIFTI_MRS.
    """
    from fsl_mrs.utils.mrs_io import read_FID
    from fsl_mrs.utils.preproc import nifti_mrs_proc

    nifti_mrs = read_FID(in_file)
    
    if dim in nifti_mrs.dim_tags:
        print("Performing align in DIM_DYN.")
        nifti_mrs = nifti_mrs_proc.align(nifti_mrs, dim)

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

def _eddy_current_correction(in_file, out_file, ecc_ref):
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
        self.inputs.out_file = _average(
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

def _phase_correct(NIFTI_MRS, ppmlim=(2.9,3.1)):
    """
    Phase correct NIFTI MRS.
    """
    from fsl_mrs.utils.preproc import NIFTI_MRS_proc
    from fsl_mrs.utils.mrs_io import read_FID
    
    NIFTI_MRS = NIFTI_MRS_proc.phase_correct(NIFTI_MRS, ppmlim=ppmlim)
    
    return NIFTI_MRS

class PhaseCorrect(BaseInterface):
    input_spec = PhaseCorrectInputSpec
    output_spec = Base_NIFTI_MRS_OutputSpec
    
    def _run_interface(self, runtime):
        # Lógica que procesa el objeto de entrada
        self._main_output = self._phase_correct_svs(
            self.inputs.main_input,
            self.inputs.ppmlim
            )
        return runtime

    def _list_outputs(self):
        # Retorna la instancia procesada
        outputs = self._outputs().get()
        outputs["main_output"] = self._main_output
        return outputs


class Shift2Creatine(BaseInterface):
    input_spec = Base_NIFTI_MRS_InputSpec
    output_spec = NIFTI_MRS_OutputSpec
    
    def _shift_to_creatine(NIFTI_MRS):
        """
        Shift to Creatine 3.027 ppm peak.
        """
        from fsl_mrs.utils.preproc import NIFTI_MRS_proc
        
        NIFTI_MRS = NIFTI_MRS_proc.shift_to_reference(
            NIFTI_MRS, 
            ppm_ref=3.027, 
            peak_search=(2.9, 3.1)
            )
        
        return NIFTI_MRS
    
    def _run_interface(self, runtime):
        # Lógica que procesa el objeto de entrada
        self._main_output = self._shift_to_creatine(self.inputs.main_input)
        return runtime

    def _list_outputs(self):
        # Retorna la instancia procesada
        outputs = self._outputs().get()
        outputs["main_output"] = self._main_output
        return outputs
    
    
class RemoveWater(BaseInterface):
    input_spec = Base_NIFTI_MRS_InputSpec
    output_spec = NIFTI_MRS_OutputSpec
    
    def _remove_water(NIFTI_MRS):
        """
        Perform eddy current correction.
        """
        from fsl_mrs.utils.preproc import NIFTI_MRS_proc
        
        hlsvdlimits = [-0.25, 0.25]
        NIFTI_MRS = NIFTI_MRS_proc.remove_peaks(NIFTI_MRS)
        
        return NIFTI_MRS
    
    def _run_interface(self, runtime):
        # Lógica que procesa el objeto de entrada
        self._main_output = self._remove_water(self.inputs.main_input)
        return runtime

    def _list_outputs(self):
        # Retorna la instancia procesada
        outputs = self._outputs().get()
        outputs["main_output"] = self._main_output
        return outputs


class DynamicAlign(BaseInterface):
    class Dynamic_NIFTI_MRS_InputSpec(
        Dim_NIFTI_MRS_InputSpec, 
        Model_NIFTI_MRS_InputSpec, 
        ppm_NIFTI_MRS_InputSpec):
        pass
    
    input_spec = Dynamic_NIFTI_MRS_InputSpec
    output_spec = NIFTI_MRS_OutputSpec
    
    def _dynamic_align(NIFTI_MRS, dim, basis, baseline_order, model, ppmlim):
        """
        Perform eddy current correction.
        """
        from fsl_mrs.utils.preproc import dyn_based_proc as dproc

        fitargs = {
            "baseline_order":baseline_order, 
            "model":model, 
            "ppmlim":ppmlim
            }

        if dim in NIFTI_MRS.dim_tags:
            NIFTI_MRS, _eps, _phi = dproc.align_by_dynamic_fit(NIFTI_MRS, basis, fitargs)
        
        return NIFTI_MRS
    
    def _run_interface(self, runtime):
        # Lógica que procesa el objeto de entrada
        self._main_output = self._dynamic_align(
            self.inputs.main_input, 
            self.inputs.dim,
            self.inputs.basis,
            self.inputs.baseline_order,
            self.inputs.model,
            self.inputs.ppmlim
            )
        return runtime

    def _list_outputs(self):
        # Retorna la instancia procesada
        outputs = self._outputs().get()
        outputs["main_output"] = self._main_output
        return outputs


# Workflows

ref_align_dyn = Node(Align(dim='DIM_DYN'), name="align ref DIM_DYN")
ref_align_edit = Node(Align(dim='DIM_EDIT'), name="align ref DIM_EDIT")
ref_average_dyn = Node(Average(dim='DIM_DYN'), name="average ECC ref dyn")
ref_average_edit = Node(Average(dim='DIM_EDIT'), name="average ECC ref edit")
ref_ecc = Node(EddyCurrentCorrection, name="ecc ref")
ref_phase_correct = Node(PhaseCorrect(ppmlim), name="phase correct ref")

ecc_ref_average_dyn = Node(Average(dim='DIM_DYN'), name="average ECC ref dyn")
ecc_ref_average_edit = Node(Average(dim='DIM_EDIT'), name="average ECC ref edit")

svs_align_dyn = Node(Align(dim='DIM_DYN'), name="align svs DIM_DYN")
svs_align_edit = Node(Align(dim='DIM_DYN'), name="align svs DIM_EDIT")
svs_average_dyn = Node(Average(dim='DIM_DYN'), name="average svs DIM_DYN")
svs_ecc = Node(EddyCurrentCorrection, name="ecc svs")
svs_phase_correct = Node(PhaseCorrect(ppmlim), name="phase correct svs")
svs_remove_water = Node(RemoveWater(), name="remove water svs")

hermes_sort = Workflow()
hermes_sort.connect(
   (svs_average_dyn, svs_remove_water, [("main_output", "main_input")]),
   (svs_remove_water, shift2creatine, [("main_output", "main_input")]),
   (shift2creatine, svs_phase_correct, [("main_output", "main_input")]),
   (svs_phase_correct, hermes_yshift, [("main_output", "main_input")]),
   (hermes_yshift, hermes_subspectra_sort, [("main_output", "main_input")]),
)

water_proc = Workflow()
water_proc.connect(
    # align water ref
    (ref_align_dyn, ref_align_edit, [("main_output", "main_input")]),
    
    # get ecc ref
    (ref_align_edit, ecc_ref_average_dyn, [("main_output", "main_input")]),
    (ecc_ref_average_dyn, ecc_ref_average_edit, [("main_output", "main_input")]),
    
    # do ecc on water ref
    (ref_align_edit, ref_ecc, [("main_output", "main_input")]),
    (ecc_ref_average_edit, ref_ecc, [("main_output", "ref")]),
    
    # phase correct water
    (ref_ecc, ref_phase_correct ,[("main_output", "main_input")]),
    
    # average water
    (ref_phase_correct, ref_average_dyn ,[("main_output", "main_input")]),
    (ref_average_dyn, ref_average_edit,[("main_output", "main_input")]),
)


hermes_proc.connect(

    # svs ecc
    (svs_align_dyn, svs_ecc, [("main_output", "main_input")]),
    (water_proc, svs_ecc, [("main_output", "ref")]),

    # sort subspectra
    (svs_ecc, hermes_sort, [("main_output", "main_input")]),
    
    # dynamic align
    (hermes_sort, svs_dynamic_align_dyn, [("main_output", "main_input")]),
    
    # remove water
    (svs_dynamic_align_dyn, svs_remove_water, [("main_output", "main_input")]),
    
    # shift to creatine peak
    (svs_remove_water, shift2creatine, [("main_output", "main_input")]),
    
    # phase correction
    (shift2creatine, svs_phase_correct, [("main_output", "main_input")]),
    
    # dim_edit y-shift
    (svs_phase_correct, hermes_yshift, [("main_output", "main_input")]),
    
    # average dim_dyn
    (hermes_yshift, svs_average_dyn, [("main_output", "main_input")]),
    
    # align dim_edit
    (svs_average_dyn, svs_align_edit, [("main_output", "main_input")]),
    
    # get edit results
    (svs_align_edit, hermes_sum, [("main_output", "main_input")]),
)