from nipype.pipeline.engine import Node, Workflow
from nipype.interfaces.base import BaseInterface, BaseInterfaceInputSpec, TraitedSpec, traits, File
from nipype.interfaces.utility import IdentityInterface, Merge
from nipype.interfaces.io import DataSink, SelectFiles
from nipype.interfaces import fsl
from mrs_pipe.interfaces.fsl_mrs import *
from mrs_pipe.interfaces.osprey import SpectralRegistration
from fsl_mrs.utils import mrs_io

def get_svs_anat_wf():

    def get_fslanat_derivatives(fsl_anat_dir):
        from pathlib import Path
        t1w2mni_file = str(Path(fsl_anat_dir) / 'T1_to_MNI_lin.nii.gz')
        t1w2mni_mat = str(Path(fsl_anat_dir) / 'T1_to_MNI_lin.mat')

        return t1w2mni_file, t1w2mni_mat
    
    Get_fslanat_derivatives = Function(
        input_names='fsl_anat_dir', 
        output_names=['t1w2mni_file', 't1w2mni_mat'], 
        function=get_fslanat_derivatives
        )

    wf = Workflow(name='svs_anat')
    inputnode = Node(IdentityInterface(fields=['t1w', 'fsl_anat_dir', 'svs']), 'inputnode')
    svs_mask = Node(svs_segment_Interface(), name='svs_mask')
    mask2standard = Node(fsl.FLIRT(apply_xfm=True, interp='nearestneighbour', uses_qform=True), 'mask2standard')
    fslanat_derivatives = Node(Get_fslanat_derivatives, name='get_fslanat_derivatives')

    wf.connect([
        (inputnode, svs_mask, [('svs', 'svs')]),
        (inputnode, svs_mask, [('fsl_anat_dir', 'fsl_anat_dir')]),
        (inputnode, fslanat_derivatives, [('fsl_anat_dir', 'fsl_anat_dir')]),
        (fslanat_derivatives, mask2standard, [('t1w2mni_file', 'reference')]),
        (fslanat_derivatives, mask2standard, [('t1w2mni_mat', 'in_matrix_file')]),
        (svs_mask, mask2standard, [('mask', 'in_file')]),
    ])
    
    return wf

def get_press_proc_wf(svs, mrsref=None, eccref=None):
    merge_svs = Node(Merge(numinputs=10), name='svs_reports2list')
    merge_mrsref = Node(Merge(numinputs=10), name='mrsref_reports2list')

    wf = Workflow(name='press_proc', base_dir='work')
    inputnode = Node(IdentityInterface(fields=['svs','mrsref', 'eccref']), name='inputnode')

    # eccref preproc
    if eccref:

        eccref_dim_tags = mrs_io.read_FID(eccref).dim_tags
        if 'DIM_DYN' in eccref_dim_tags:

            # frequency and phase align eccref transients (DIM_DYN)
            align_eccref = Node(Align(), name='align_eccref')
            align_eccref.inputs.dim = 'DIM_DYN'
            align_eccref.inputs.ppmlim = (0,8)
        
            # average 
            average_eccref = Node(Average(), name='avg_dyn_eccref')
            average_eccref.inputs.dim = 'DIM_DYN'

            wf.connect([
                (inputnode, align_eccref, [('mrsref', 'in_file')]),
                (align_eccref, average_eccref, [('out_file', 'in_file')])
            ])

    # mrsref preproc
    if mrsref:

        mrsref_dim_tags = mrs_io.read_FID(mrsref).dim_tags

        mrsref_merge_mrs_reports = Node(merge_mrs_reports_Interface(), name='mrsref_merge_mrs_reports')
        get_mrsref_stem = Node(get_file_stem_Interface, 'get_mrsref_stem')

        # phase correction
        phase_correct_ref = Node(PhaseCorrect(), name='phase_correct_mrsref')
        phase_correct_ref.inputs.ppmlim = (4.55, 4.7)
        phase_correct_ref.inputs.out_file = 'preproc'

        ecc_mrsref = Node(EddyCurrentCorrection(), name='ec_correct_mrsref')
        
        if 'DIM_DYN' in mrsref_dim_tags:

            # frequency and phase align mrsref transients (DIM_DYN)
            align_mrsref = Node(Align(), name='align_mrsref')
            align_mrsref.inputs.dim = 'DIM_DYN'
            align_mrsref.inputs.ppmlim = (0,8)
        
            # average 
            average_mrsref = Node(Average(), name='avg_dyn_mrsref')
            average_mrsref.inputs.dim = 'DIM_DYN'

            wf.connect([
                (inputnode, align_mrsref, [('mrsref', 'in_file')]),
                (align_mrsref, average_mrsref, [('out_file', 'in_file')]),
                (average_mrsref, ecc_mrsref, [('out_file', 'in_file')])

                (align_mrsref, merge_mrsref, [('report', 'in1')]),
                (average_mrsref, merge_mrsref, [('report', 'in2')])
            ])

        else:
            wf.connect([
                (inputnode, ecc_mrsref, [('mrsref', 'in_file')])
            ])

        # perform ecc on mrsref
        if eccref:
            if 'DIM_DYN' in eccref_dim_tags:
                wf.connect([
                    (average_eccref, ecc_mrsref, [('out_file', 'ref')])
                ])
            else:
                wf.connect([
                    (inputnode, ecc_mrsref, [('eccref', 'ref')])
                ])
        else:
            if 'DIM_DYN' in mrsref_dim_tags:
                wf.connect([
                    (average_mrsref, ecc_mrsref, [('out_file', 'ref')])
                ])
            else:
                wf.connect([
                    (inputnode, ecc_mrsref, [('mrsref', 'ref')])
                ])
    
        # perform phase correction
        wf.connect([
            (ecc_mrsref, phase_correct_ref, [('out_file', 'in_file')]),

            (ecc_mrsref, merge_mrsref, [('report', 'in3')]),
            (phase_correct_ref, merge_mrsref, [('report', 'inf4')]),
            
            (merge_mrsref, mrsref_merge_mrs_reports, [('out', 'in_files')]),
            (inputnode, get_mrsref_stem, [('mrsref', 'file')]),
            (get_mrsref_stem, mrsref_merge_mrs_reports, [('out', 'out_file')]),
            (get_mrsref_stem, mrsref_merge_mrs_reports, [('out', 'description')])
        ])


    # svs preprocessing

    remove_water = Node(RemoveWater(), name='remove_water')
    phase_correct_svs = Node(PhaseCorrect_Creatine_ppmlim(), name='phase_correct_svs')
    shift2creatine = Node(ShiftToCreatine(), name='shift2creatine')
    svs_dim_tags = mrs_io.read_FID(svs).dim_tags

    ## align transients
    if 'DIM_DYN' in svs_dim_tags:

        align_dyn = Node(Align(), name='align_dyn')
        align_dyn.inputs.dim = 'DIM_DYN'

        wf.connect([
            (inputnode, align_dyn, [('svs', 'in_file')]),
            (align_dyn, merge_svs, [('report', 'in1')])
        ])

    ## ecc svs
    if inputnode.inputs.eccref or inputnode.inputs.mrsref:
        ecc_svs = Node(EddyCurrentCorrection(), name='ecc_svs')

        if 'DIM_DYN' in svs_dim_tags:
            wf.connect([
                (align_dyn, ecc_svs, [('out_file', 'in_file')])
            ])
        else:
            wf.connect([
                (inputnode, ecc_svs, [('svs', 'in_file')])
            ])

        if inputnode.inputs.eccref:
            if 'DIM_DYN' in eccref_dim_tags:
                wf.connect([
                    (average_eccref, ecc_svs, [('out_file', 'ref')])
                ])
            else:
                wf.connect([
                    (inputnode, ecc_svs, [('eccref', 'ref')])
                ])
        elif inputnode.inputs.mrsref:
            if 'DIM_DYN' in mrsref_dim_tags:
                wf.connect([
                    (average_mrsref, ecc_svs, [('out_file', 'ref')])
                ])
            else:
                wf.connect([
                    (inputnode, ecc_svs, [('mrsref', 'ref')])
                ])

    ## remove water
        wf.connect([
            (ecc_svs, remove_water, [('out_file', 'in_file')]),
            (ecc_svs, merge_svs, [('report', 'in2')])
        ])

    else:
        if 'DIM_DYN' in svs_dim_tags:
            wf.connect([
                (align_dyn, remove_water, [('out_file', 'in_file')])
            ])
        else:
            wf.connect([
                (inputnode, remove_water, [('svs', 'in_file')])
            ])

    wf.connect([
        (remove_water, shift2creatine, [('out_file', 'in_file')]),
        (shift2creatine, phase_correct_svs,  [('out_file', 'in_file')]),

        (remove_water, merge_svs, [('report', 'in3')]),
        (shift2creatine, merge_svs, [('report', 'in4')]),
        (phase_correct_svs, merge_svs, [('report', 'in5')])
    ])

    if 'DIM_DYN' in svs_dim_tags:
        average_svs = Node(Average(), name='avg_dyn')
        average_svs.inputs.dim = 'DIM_DYN'
        average_svs.inputs.out_file = 'preproc'
    
        wf.connect([
            (phase_correct_svs, average_svs, [('out_file', 'in_file')]),
            (average_svs, merge_svs, [('report', 'in6')])
        ])
    else:
        phase_correct_svs.inputs.out_file = 'preproc'

    svs_merge_mrs_reports = Node(merge_mrs_reports_Interface(), name='svs_merge_mrs_reports')
    get_svs_stem = Node(get_file_stem_Interface, 'get_svs_stem')
    
    wf.connect([
        (inputnode, get_svs_stem,  [('svs', 'file')]),
        (get_svs_stem, svs_merge_mrs_reports, [('out', 'description')]),
        (get_svs_stem, svs_merge_mrs_reports, [('out', 'out_file')]),
        (merge_svs, svs_merge_mrs_reports, [('out', 'in_files')])
    ])
    
    
    return wf

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
    # average_svs = Node(Average(), name='avg_dyn')
    # average_svs.inputs.dim = 'DIM_DYN'
    align_edit = Node(Align(), name='align_edit')
    align_edit.inputs.dim = 'DIM_EDIT'
    align_edit.inputs.out_file = 'preproc'
    align_edit.inputs.report_append = 'DIM_EDIT'
    hermes_yalign = Node(HERMESAlignYEdit(), name='hermes_yalign')
    prob_specreg = Node(SpectralRegistration(), name='prob_specreg')
    prob_specreg.inputs.osprey_path='/misc/geminis2/ramirezd/osprey/'
    prob_specreg.inputs.method='Probabilistic'
    prob_specreg.inputs.sequence='HERMES'
    prob_specreg.inputs.spm_path='/misc/geminis2/spm12/'
    prob_specreg.inputs.out_file='preproc'
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
        (align_edit, prob_specreg, [('out_file', 'in_file')]),
        (prob_specreg, edit_sum, [('out_file', 'in_file')]),

        (infosource, preproc_sinker,  [('subject', 'container')]),

        (edit_sum, preproc_sinker, [('summation', 'mrs.@summation')]),
        (edit_sum, preproc_sinker, [('gaba', 'mrs.@gaba')]),
        (edit_sum, preproc_sinker, [('gsh', 'mrs.@gsh')]),

        (prob_specreg, preproc_sinker, [('out_file', 'mrs.@average_svs')]),
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


