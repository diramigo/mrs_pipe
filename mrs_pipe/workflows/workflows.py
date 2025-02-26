from nipype.pipeline.engine import Node, Workflow
from nipype.interfaces.base import BaseInterface, BaseInterfaceInputSpec, TraitedSpec, traits, File
from nipype.interfaces.utility import IdentityInterface, Merge
from nipype.interfaces.io import DataSink, SelectFiles
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

def get_hermes_subspectra_proc_wf(name, report_append=None):

    align_dyn = Node(Align(), name='align_dyn')
    align_dyn.inputs.dim = 'DIM_DYN'
    ec_correct_svs = Node(EddyCurrentCorrection(), name='ec_correct_svs')
    remove_water = Node(RemoveWater(), name='remove_water')
    shift2creatine = Node(ShiftToCreatine(), name='shift2creatine')
    phase_correct_svs = Node(PhaseCorrect_Creatine_ppmlim(), name='phase_correct_svs')
    
    if report_append:
        align_dyn.inputs.report_append = report_append
        ec_correct_svs.inputs.report_append = report_append
        remove_water.inputs.report_append = report_append
        shift2creatine.inputs.report_append = report_append
        phase_correct_svs.inputs.report_append = report_append

    svs_reports2list = Node(Merge(numinputs=5), name='svs_reports2list')
    hermes_subspectra_proc_wf = Workflow(name=name, base_dir='work')

    hermes_subspectra_proc_wf.connect([
        (align_dyn, ec_correct_svs, [('out_file', 'in_file')]),
        (ec_correct_svs, remove_water, [('out_file', 'in_file')]),
        (remove_water, shift2creatine, [('out_file', 'in_file')]),
        (shift2creatine, phase_correct_svs, [('out_file', 'in_file')]),

        (align_dyn, svs_reports2list, [('report', 'in1')]),
        (ec_correct_svs, svs_reports2list, [('report', 'in2')]),
        (remove_water, svs_reports2list, [('report', 'in3')]),
        (shift2creatine, svs_reports2list, [('report', 'in4')]),
        (phase_correct_svs, svs_reports2list, [('report', 'in5')]),
    ])

    return hermes_subspectra_proc_wf

def get_hermes_proc_wf(subject, voi, bids_root, output_dir):

    bids_root = os.path.abspath(bids_root)
    output_dir = os.path.abspath(output_dir)

    infosource = Node(IdentityInterface(fields=["subject", "voi"]),
                    name="infosource")
    infosource.iterables = [("subject", subject), ("voi", voi)]

    # select infiles
    templates = {
        "svs": "{subject}/mrs/{subject}_acq-hermes_voi-{voi}_svs.nii.gz",
        "ref": "{subject}/mrs/{subject}_acq-hermes_voi-{voi}_mrsref.nii.gz"
    }
    select_files = Node(SelectFiles(templates), 'select_files')
    select_files.inputs.base_directory = bids_root

    # svs
    get_svs_stem = Node(get_file_stem_Interface, 'get_svs_stem')

    hermes_sort_subspectra = Node(HERMESSortSubspectra(), name='hermes_sort_subspectra')

    split_hermes = Node(SplitHermes(), name='split_hermes')
    proc_gaba_on = get_hermes_subspectra_proc_wf('proc_gaba_on', 'GABA_ON')
    proc_gsh_on = get_hermes_subspectra_proc_wf('proc_gsh_on', 'GSH_ON')
    proc_both_off = get_hermes_subspectra_proc_wf('proc_both_off', 'BOTH_OFF')
    proc_both_on = get_hermes_subspectra_proc_wf('proc_both_on', 'BOTH_ON')

    merge_hermes = Node(MergeHermes(), name='merge_hermes')
    average_svs = Node(Average(), name='avg_dyn')
    average_svs.inputs.dim = 'DIM_DYN'
    align_edit = Node(Align(), name='align_edit')
    align_edit.inputs.dim = 'DIM_EDIT'
    align_edit.inputs.out_file = 'preproc'
    align_edit.inputs.report_append = 'DIM_EDIT'
    hermes_yalign = Node(HERMESAlignYEdit(), name='hermes_yalign')
    edit_sum = Node(HERMES_Edit_Sum(), name='edit_sum')

    ## merge svs reports 
    svs_reports2list = Node(Merge(numinputs=5), name='svs_reports2list')
    svs_merge_mrs_reports = Node(merge_mrs_reports_Interface(), name='svs_merge_mrs_reports')

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

    ref_reports2list = Node(Merge(numinputs=4), name='ref_reports2list')
    ref_merge_mrs_reports = Node(merge_mrs_reports_Interface(), name='ref_merge_mrs_reports')
    get_ref_stem = Node(get_file_stem_Interface, 'get_ref_stem')

    # sink
    preproc_sinker = Node(DataSink(), name='preproc_derivatives')
    preproc_sinker.inputs.base_directory = output_dir
    preproc_sinker.inputs.parameterization = False

    report_sinker = Node(
        DataSink(
            base_directory=output_dir, 
            container='', 
            parameterization=False), 
        name='report_derivatives'
        )


    hermes_proc = Workflow(name='hermes_proc', base_dir='work')

    hermes_proc.connect(
        [
        (infosource, select_files,  [('subject', 'subject')]),
        (infosource, select_files,  [('voi', 'voi')]),
        # ref
        (select_files, get_ref_stem,  [('ref', 'file')]),
        (select_files, remove_zero_mean_ref,  [('ref', 'in_file')]),
        (remove_zero_mean_ref, align_ref,  [('out_file', 'in_file')]),
        (align_ref, get_ecc_ref,  [('out_file', 'in_file')]),
        (get_ecc_ref, ec_correct_ref,  [('out_file', 'ref')]),
        (align_ref, ec_correct_ref,  [('out_file', 'in_file')]),
        (ec_correct_ref, phase_correct_ref, [('out_file', 'in_file')]),
        (phase_correct_ref, average_ref, [('out_file', 'in_file')]),

        #svs
        (select_files, hermes_sort_subspectra,  [('svs', 'in_file')]),
        (select_files, get_svs_stem,  [('svs', 'file')]),
        (hermes_sort_subspectra, split_hermes,  [('out_file', 'in_file')]),
        (split_hermes, proc_gaba_on,  [('gaba_on', 'align_dyn.in_file')]),
        (get_ecc_ref, proc_gaba_on, [('out_file', 'ec_correct_svs.ref')]),
        (split_hermes, proc_gsh_on,  [('gsh_on', 'align_dyn.in_file')]),
        (get_ecc_ref, proc_gsh_on, [('out_file', 'ec_correct_svs.ref')]),
        (split_hermes, proc_both_on,  [('both_on', 'align_dyn.in_file')]),
        (get_ecc_ref, proc_both_on, [('out_file', 'ec_correct_svs.ref')]),
        (split_hermes, proc_both_off,  [('both_off', 'align_dyn.in_file')]),
        (get_ecc_ref, proc_both_off, [('out_file', 'ec_correct_svs.ref')]),

        (proc_gaba_on, merge_hermes,  [('phase_correct_svs.out_file', 'gaba_on')]),
        (proc_gsh_on, merge_hermes,  [('phase_correct_svs.out_file', 'gsh_on')]),
        (proc_both_on, merge_hermes,  [('phase_correct_svs.out_file', 'both_on')]),
        (proc_both_off, merge_hermes,  [('phase_correct_svs.out_file', 'both_off')]),

        (merge_hermes, hermes_yalign, [('out_file', 'in_file')]),
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
        (phase_correct_ref, preproc_sinker, [('out_file', 'mrs.@preproc_ref')]),

        ## reports svs
        (proc_gaba_on, svs_reports2list, [('svs_reports2list.out', 'in1')]),
        (proc_gsh_on, svs_reports2list, [('svs_reports2list.out', 'in2')]),
        (proc_both_on, svs_reports2list, [('svs_reports2list.out', 'in3')]),
        (proc_both_off, svs_reports2list, [('svs_reports2list.out', 'in4')]),
        (align_edit, svs_reports2list, [('report', 'in5')]),
        (get_svs_stem, svs_merge_mrs_reports, [('out', 'out_file')]),
        (get_svs_stem, svs_merge_mrs_reports, [('out', 'description')]),
        (svs_reports2list, svs_merge_mrs_reports, [('out', 'in_files')]),
        (svs_merge_mrs_reports, report_sinker, [('out_file', '@svs_report')]),

        ## reports ref
        (align_ref, ref_reports2list,  [('report', 'in1')]),
        (average_ref, ref_reports2list,  [('report', 'in2')]),
        (ec_correct_ref, ref_reports2list, [('report', 'in3')]),
        (phase_correct_ref, ref_reports2list, [('report', 'in4')]),
        (get_ref_stem, ref_merge_mrs_reports, [('out', 'out_file')]),
        (get_ref_stem, ref_merge_mrs_reports, [('out', 'description')]),
        (ref_reports2list, ref_merge_mrs_reports, [('out', 'in_files')]),
        (ref_merge_mrs_reports, report_sinker, [('out_file', '@ref_report')]),

        ]

    )

    return hermes_proc


