from nipype.interfaces.base import BaseInterface, BaseInterfaceInputSpec, OutputSpec, traits
from fsl_mrs.utils.nifti_mrs import NIfTI_MRS
from fsl_mrs.utils.mrs_basis import Basis

# Functions

def read_nifti_mrs(path):
    """
    Read NIfTI MRS
    """
    
    from fsl_mrs.utils import mrs_io
    
    nifti_mrs = mrs_io.read_FID(path)
    
    return nifti_mrs


def read_basis(path, type=['press', 'hermes']):
    """
    Read NIfTI MRS
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
            c = mrs_io.basis(Path(path) / 'd')
            return a,b,c,d
        case _:
            return None
        

def preproc_spectrum(svs, ecc_ref=None):
    """
    Preprocess HERMES spectra that have been aligned within DIM_DYN.
    """
    if ecc_ref is not None:
        svs = nifti_mrs_proc.ecc(svs, ecc_ref)

    hlsvdlimits = [-0.25, 0.25]
    svs = nifti_mrs_proc.remove_peaks(svs, hlsvdlimits, limit_units='ppm')
    svs = nifti_mrs_proc.shift_to_reference(svs, ppm_ref=3.027, peak_search=(2.9, 3.1))
    svs = nifti_mrs_proc.phase_correct(svs, ppmlim=(2.9, 3.1))
    return svs


# Input and Output Specs

class Base_NIfTI_MRS_InputSpec(BaseInterfaceInputSpec):
    main_input = traits.Instance(
        NIfTI_MRS, 
        desc="Main NIfTI_MRS input", 
        mandatory=True
        )
    
class Ref_NIfTI_MRS_InputSpec(Base_NIfTI_MRS_InputSpec):
    ref = traits.Instance(
        NIfTI_MRS,
        desc="NIfTI_MRS water reference input", 
        mandatory=True
        )

class Dim_NIfTI_MRS_InputSpec(Base_NIfTI_MRS_InputSpec):
    dim = traits.Enum(
        "DIM_DYN", "DIM_EDIT",
        desc="NIfTI_MRS dimension ('DIM_DYN':transients, 'DIM_EDIT':edit subspectra).",
        mandatory=True
    )
    
class ppm_NIfTI_MRS_InputSpec(Base_NIfTI_MRS_InputSpec):
    ppmlim = traits.Tuple(traits.Float, desc="ppm range.", mandatory=True)
    
class Model_NIfTI_MRS_InputSpec(Base_NIfTI_MRS_InputSpec):
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
    
    # num_threads = traits.Int(
    #     "Numbers of threads when running interface.", 
    #     default=1, 
    #     usedefault=False
    #     )

class NIfTI_MRS_OutputSpec(OutputSpec):
    main_output = traits.Either(
        traits.Instance(NIfTI_MRS),
        traits.Tuple(traits.Instance(NIfTI_MRS)),
        desc="Main NIfTI_MRS output"
        )

# Interfaces

class Align(BaseInterface):
    input_spec = Dim_NIfTI_MRS_InputSpec
    output_spec = NIfTI_MRS_OutputSpec
    
    def _align(nifti_mrs, dim='DIM_DYN'):
        """
        Align NIfTI_MRS.
        """
        from fsl_mrs.utils.preproc import nifti_mrs_proc
        
        if dim in nifti_mrs.dim_tags:
            nifti_mrs = nifti_mrs_proc.align(nifti_mrs, dim)
        
        return nifti_mrs
    
    def _run_interface(self, runtime):
        # Lógica que procesa el objeto de entrada
        self._main_output = self._align(
            self.inputs.main_input, 
            self.inputs.dim
            )
        return runtime

    def _list_outputs(self):
        # Retorna la instancia procesada
        outputs = self._outputs().get()
        outputs["main_output"] = self._main_output
        return outputs
    
    
class Average(BaseInterface):
    input_spec = Dim_NIfTI_MRS_InputSpec
    output_spec = NIfTI_MRS_OutputSpec
    
    def _average(nifti_mrs, dim='DIM_DYN'):
        """
        Average NIfTI_MRS.
        """
        from fsl_mrs.utils.preproc import nifti_mrs_proc
        
        if dim in ref.dim_tags:
            nifti_mrs = nifti_mrs_proc._average(nifti_mrs, dim)
            
        return nifti_mrs
    
    def _run_interface(self, runtime):
        # Lógica que procesa el objeto de entrada
        self._main_output = self._average(
            self.inputs.main_input,
            self.inputs.dim
            )
        return runtime

    def _list_outputs(self):
        # Retorna la instancia procesada
        outputs = self._outputs().get()
        outputs["main_output"] = self._main_output
        return outputs


class EddyCurrentCorrection(BaseInterface):
    input_spec = Ref_NIfTI_MRS_InputSpec
    output_spec = NIfTI_MRS_OutputSpec
    
    def _eddy_current_correction(nifti_mrs, ecc_ref=None):
        """
        Perform eddy current correction.
        """
        from fsl_mrs.utils.preproc import nifti_mrs_proc
        
        if ecc_ref is not None:
            nifti_mrs = nifti_mrs_proc.ecc(nifti_mrs, ecc_ref)
        
        return nifti_mrs
    
    def _run_interface(self, runtime):
        # Lógica que procesa el objeto de entrada
        self._main_output = self._eddy_current_correction(
            self.inputs.main_input,
            self.inputs.ref
            )
        return runtime

    def _list_outputs(self):
        # Retorna la instancia procesada
        outputs = self._outputs().get()
        outputs["main_output"] = self._main_output
        return outputs


class PhaseCorrect(BaseInterface):
    input_spec = ppm_NIfTI_MRS_InputSpec
    output_spec = NIfTI_MRS_OutputSpec
    
    def _phase_correct(nifti_mrs, ppmlim=(2.9,3.1)):
        """
        Phase correct NIfTI MRS.
        """
        from fsl_mrs.utils.preproc import nifti_mrs_proc
        
        nifti_mrs = nifti_mrs_proc.phase_correct(nifti_mrs, ppmlim=ppmlim)
        
        return nifti_mrs
    
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
    input_spec = Base_NIfTI_MRS_InputSpec
    output_spec = NIfTI_MRS_OutputSpec
    
    def _shift_to_creatine(nifti_mrs):
        """
        Shift to Creatine 3.027 ppm peak.
        """
        from fsl_mrs.utils.preproc import nifti_mrs_proc
        
        nifti_mrs = nifti_mrs_proc.shift_to_reference(
            nifti_mrs, 
            ppm_ref=3.027, 
            peak_search=(2.9, 3.1)
            )
        
        return nifti_mrs
    
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
    input_spec = Base_NIfTI_MRS_InputSpec
    output_spec = NIfTI_MRS_OutputSpec
    
    def _remove_water(nifti_mrs):
        """
        Perform eddy current correction.
        """
        from fsl_mrs.utils.preproc import nifti_mrs_proc
        
        hlsvdlimits = [-0.25, 0.25]
        nifti_mrs = nifti_mrs_proc.remove_peaks(nifti_mrs)
        
        return nifti_mrs
    
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
    class Dynamic_NIfTI_MRS_InputSpec(
        Dim_NIfTI_MRS_InputSpec, 
        Model_NIfTI_MRS_InputSpec, 
        ppm_NIfTI_MRS_InputSpec):
        pass
    
    input_spec = Dynamic_NIfTI_MRS_InputSpec
    output_spec = NIfTI_MRS_OutputSpec
    
    def _dynamic_align(nifti_mrs, dim, basis, baseline_order, model, ppmlim):
        """
        Perform eddy current correction.
        """
        from fsl_mrs.utils.preproc import dyn_based_proc as dproc

        fitargs = {
            "baseline_order":baseline_order, 
            "model":model, 
            "ppmlim":ppmlim
            }

        if dim in nifti_mrs.dim_tags:
            nifti_mrs, _eps, _phi = dproc.align_by_dynamic_fit(nifti_mrs, basis, fitargs)
        
        return nifti_mrs
    
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