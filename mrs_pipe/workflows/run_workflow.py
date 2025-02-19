if __name__ == '__main__':

    import os
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-s','--subject_label', nargs='+')
    parser.add_argument('-v','--volume', nargs='+')
    parser.add_argument('--n_procs', type=int, default=1)
    parser.add_argument('--omp_nthreads', type=int, default=1)
    parser.add_argument('bids_dir', type=str)
    parser.add_argument('output_dir', type=str)
    args = parser.parse_args()

    # set omp_nthreads before importing packages with numpy dependencies
    if args.omp_nthreads:
        os.environ["OMP_NUM_THREADS"] = str(args.omp_nthreads)

    from mrs_pipe.workflows.workflows import get_hermes_proc_wf

    basic_proc = get_hermes_proc_wf(args.subject_label, args.volume, args.bids_dir, args.output_dir)
    basic_proc.run(plugin='MultiProc',  plugin_args={'n_procs' : args.n_procs})