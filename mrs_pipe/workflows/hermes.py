from nipype.pipeline.engine import Node, Workflow
from nipype.interfaces.base import BaseInterface, BaseInterfaceInputSpec, TraitedSpec, traits, File

from mrs_pipe.interfaces import fsl_mrs 


def get_basic_proc_wf():

    # svs
    align_dyn = Node(Align(), name='align_dyn')
    align_dyn.inputs.dim = 'DIM_DYN'
    align_dyn.inputs.out_file = 'svs_aligned.nii.gz'

    align_edit = Node(Align(), name='align_edit')
    align_edit.inputs.dim = 'DIM_EDIT'
    align_edit.inputs.out_file = 'svs_aligned.nii.gz'

    average = Node(Average(), name='avg_dyn')
    average.inputs.dim = 'DIM_DYN'
    average.inputs.out_file = 'svs_averaged.nii.gz'

    ec_correct_svs = Node(EddyCurrentCorrection(), name='ec_correct_svs')
    ec_correct_svs.inputs.out_file = 'svs_ec_corrected.nii.gz'

    remove_water = Node(RemoveWater(), name='remove_water')
    remove_water.inputs.out_file = 'svs_water_removed.nii.gz'

    shift2creatine = Node(ShiftToCreatine(), name='shift2creatine')
    shift2creatine.inputs.out_file = 'svs_shifted_cr.nii.gz'

    phase_correct_svs = Node(PhaseCorrect_Creatine_ppmlim(), name='phase_correct_svs')
    phase_correct_svs.inputs.out_file = 'svs_phase_corrected.nii.gz'

    hermes_sort_subspectra = Node(HERMESSortSubspectra(), name='hermes_sort_subspectra')
    hermes_sort_subspectra.inputs.out_file = 'svs_edit_sorted.nii.gz'
    hermes_sort_subspectra.inputs.in_file = '/home/diego/Documents/projects/mrs_gaba/inputs/mri-raw/sub-001A/mrs/sub-001A_acq-hermes_voi-pcc_svs.nii.gz'

    hermes_yalign = Node(HERMESAlignYEdit(), name='hermes_yalign')
    hermes_yalign.inputs.out_file = 'svs_yalign.nii.gz'  

    edit_sum = Node(HERMES_Edit_Sum(), name='edit_sum')

    # ref
    remove_zero_mean_ref = Node(_HERMESRefRemoveZeroMeanTransients(), name='remove_zero_mean_ref')
    remove_zero_mean_ref.inputs.in_file = '/home/diego/Documents/projects/mrs_gaba/inputs/mri-raw/sub-001A/mrs/sub-001A_acq-hermes_voi-pcc_ref.nii.gz'
    remove_zero_mean_ref.inputs.out_file = 'ref.nii.gz'

    align_ref = Node(Align(), name='align_ref')
    align_ref.inputs.dim = 'DIM_DYN'
    align_ref.inputs.ppmlim = (0,8)
    align_ref.inputs.out_file = 'ref_aligned.nii.gz'

    get_ecc_ref = Node(Average(), name='get_ecc_ref')
    get_ecc_ref.inputs.dim = 'DIM_DYN'
    get_ecc_ref.inputs.out_file = 'ecc_ref.nii.gz'

    ec_correct_ref = Node(EddyCurrentCorrection(), name='ec_correct_ref')
    ec_correct_ref.inputs.out_file = 'ref_ec_corrected.nii.gz'

    phase_correct_ref = Node(PhaseCorrect(), name='phase_correct_ref')
    phase_correct_ref.inputs.ppmlim = (4.55, 4.7)
    phase_correct_ref.inputs.out_file = 'ref_phase_corrected.nii.gz'



    basic_proc = Workflow(name='basic_proc')

    basic_proc.connect(
        [
        # ref
        (remove_zero_mean_ref, align_ref,  [('out_file', 'in_file')]),
        (align_ref, get_ecc_ref,  [('out_file', 'in_file')]),
        (get_ecc_ref, ec_correct_ref,  [('out_file', 'ref')]),
        (align_ref, ec_correct_ref,  [('out_file', 'in_file')]),
        (ec_correct_ref, phase_correct_ref, [('out_file', 'in_file')]),

        #svs
        (hermes_sort_subspectra, align_dyn,  [('out_file', 'in_file')]),
        ## ecc
        (align_dyn, ec_correct_svs, [('out_file', 'in_file')]),
        (get_ecc_ref, ec_correct_svs, [('out_file', 'ref')]),

        (ec_correct_svs, remove_water, [('out_file', 'in_file')]),
        (remove_water, shift2creatine, [('out_file', 'in_file')]),
        (shift2creatine, phase_correct_svs, [('out_file', 'in_file')]),
        (phase_correct_svs, hermes_yalign, [('out_file', 'in_file')]),
        (hermes_yalign, align_edit, [('out_file', 'in_file')]),
        (align_edit, average, [('out_file', 'in_file')]),
        (average, edit_sum, [('out_file', 'in_file')]),
        ]
    )

    return basic_proc

basic_proc = get_basic_proc_wf()
basic_proc.run()