from nipype.pipeline.engine import Node, Workflow
from nipype.interfaces.base import BaseInterface, BaseInterfaceInputSpec, TraitedSpec, traits, File
from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces.io import DataSink, SelectFiles
from nipype.interfaces import Function
from mrs_pipe.interfaces.fsl_mrs import *

def get_press_proc_wf():
    
    # svs
    phase_correct_svs = Node(
        PhaseCorrect_Creatine_ppmlim(), 
        name='phase_correct_svs'
        )
    
    align_dyn = Node(Align(), name='align_dyn')
    align_dyn.inputs.dim = 'DIM_DYN'
    
    remove_water = Node(RemoveWater(), name='remove_water')
    
    shift2creatine = Node(ShiftToCreatine(), name='shift2creatine')

    average_svs = Node(Average(), name='avg_dyn')
    average_svs.inputs.dim = 'DIM_DYN'
    
    ec_correct_svs = Node(EddyCurrentCorrection(), name='ec_correct_svs')

    # ref
    align_ref = Node(Align(), name='align_ref')
    align_ref.inputs.dim = 'DIM_DYN'
    align_ref.inputs.ppmlim = (0,8)
    
    get_ecc_ref = Node(Average(), name='get_ecc_ref')
    get_ecc_ref.inputs.dim = 'DIM_DYN'
    
    ec_correct_ref = Node(EddyCurrentCorrection(), name='ec_correct_ref')
    
    phase_correct_ref = Node(PhaseCorrect(), name='phase_correct_ref')
    phase_correct_ref.inputs.ppmlim = (4.55, 4.7)
    phase_correct_ref.inputs.out_file = 'preproc'
    
    average_ref = Node(Average(), name='avg_dyn_ref')
    average_ref.inputs.dim = 'DIM_DYN'
    
    return None


def get_hermes_proc_wf(subject, bids_root):

    infosource = Node(IdentityInterface(fields=["subject"]),
                    name="infosource")
    infosource.iterables = [("subject", subject)]

    # select infiles
    templates = {
        "svs": "{subject}/mrs/{subject}_acq-hermes_voi-pcc_svs.nii.gz",
        "ref": "{subject}/mrs/{subject}_acq-hermes_voi-pcc_ref.nii.gz"
        # "anat": "{subject}/anat/{subject}_T1w.nii.gz"

    }
    select_files = Node(SelectFiles(templates), 'select_files')
    select_files.inputs.base_directory = bids_root
    # select_files.inputs.subject = subject

    # svs
    hermes_sort_subspectra = Node(HERMESSortSubspectra(), name='hermes_sort_subspectra')
    
    align_edit = Node(Align(), name='align_edit')
    align_edit.inputs.dim = 'DIM_EDIT'
    align_edit.inputs.out_file = 'preproc'

    average_svs = Node(Average(), name='avg_dyn')
    average_svs.inputs.dim = 'DIM_DYN'
    align_dyn = Node(Align(), name='align_dyn')
    align_dyn.inputs.dim = 'DIM_DYN'
    ec_correct_svs = Node(EddyCurrentCorrection(), name='ec_correct_svs')
    remove_water = Node(RemoveWater(), name='remove_water')
    shift2creatine = Node(ShiftToCreatine(), name='shift2creatine')
    phase_correct_svs = Node(PhaseCorrect_Creatine_ppmlim(), name='phase_correct_svs')
    
    hermes_yalign = Node(HERMESAlignYEdit(), name='hermes_yalign')
    edit_sum = Node(HERMES_Edit_Sum(), name='edit_sum')

    # ref
    remove_zero_mean_ref = Node(HERMESRefRemoveZeroMeanTransients(), name='remove_zero_mean_ref')
    
    align_ref = Node(Align(), name='align_ref')
    align_ref.inputs.dim = 'DIM_DYN'
    align_ref.inputs.ppmlim = (0,8)
    
    get_ecc_ref = Node(Average(), name='get_ecc_ref')
    get_ecc_ref.inputs.dim = 'DIM_DYN'
    
    ec_correct_ref = Node(EddyCurrentCorrection(), name='ec_correct_ref')
    
    phase_correct_ref = Node(PhaseCorrect(), name='phase_correct_ref')
    phase_correct_ref.inputs.ppmlim = (4.55, 4.7)
    phase_correct_ref.inputs.out_file = 'preproc'
    
    average_ref = Node(Average(), name='avg_dyn_ref')
    average_ref.inputs.dim = 'DIM_DYN'

    # sink
    preproc_sinker = Node(DataSink(), name='preproc_derivatives')
    preproc_sinker.inputs.base_directory = os.path.abspath('./outputs')
    preproc_sinker.inputs.parameterization = False


    hermes_proc = Workflow(name='hermes_proc', base_dir='work')

    hermes_proc.connect(
        [
        (infosource, select_files,  [('subject', 'subject')]),
        # ref
        (select_files, remove_zero_mean_ref,  [('ref', 'in_file')]),
        (remove_zero_mean_ref, align_ref,  [('out_file', 'in_file')]),
        (align_ref, get_ecc_ref,  [('out_file', 'in_file')]),
        (get_ecc_ref, ec_correct_ref,  [('out_file', 'ref')]),
        (align_ref, ec_correct_ref,  [('out_file', 'in_file')]),
        (ec_correct_ref, phase_correct_ref, [('out_file', 'in_file')]),
        (phase_correct_ref, average_ref, [('out_file', 'in_file')]),

        #svs
        (select_files, hermes_sort_subspectra,  [('svs', 'in_file')]),
        (hermes_sort_subspectra, align_dyn,  [('out_file', 'in_file')]),

        (align_dyn, ec_correct_svs, [('out_file', 'in_file')]),
        (get_ecc_ref, ec_correct_svs, [('out_file', 'ref')]),

        (ec_correct_svs, remove_water, [('out_file', 'in_file')]),
        (remove_water, shift2creatine, [('out_file', 'in_file')]),
        (shift2creatine, phase_correct_svs, [('out_file', 'in_file')]),
        (phase_correct_svs, hermes_yalign, [('out_file', 'in_file')]),
        (hermes_yalign, align_edit, [('out_file', 'in_file')]),
        (align_edit, average_svs, [('out_file', 'in_file')]),
        (average_svs, edit_sum, [('out_file', 'in_file')]),

        (infosource, preproc_sinker,  [('subject', 'container')]),

        (edit_sum, preproc_sinker, [('summation', 'mrs.@summation')]),
        (edit_sum, preproc_sinker, [('gaba', 'mrs.@gaba')]),
        (edit_sum, preproc_sinker, [('gsh', 'mrs.@gsh')]),

        (average_svs, preproc_sinker, [('out_file', 'mrs.@average_svs')]),
        (average_ref, preproc_sinker, [('out_file', 'mrs.@average_ref')]),

        (align_edit, preproc_sinker, [('out_file', 'mrs.@preproc_svs')]),
        (phase_correct_ref, preproc_sinker, [('out_file', 'mrs.@preproc_ref')])

        ]
    )

    return hermes_proc


