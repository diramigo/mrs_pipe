if __name__ == '__main__':

    import os
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-s','--subject', nargs='+')
    parser.add_argument('-d','--bids_dir', type=str)
    parser.add_argument('--omp_nthreads', type=int)
    args = parser.parse_args()

    # set omp_nthreads before importing packages with numpy dependencies
    if args.omp_nthreads:
        os.environ["OMP_NUM_THREADS"] = str(args.omp_nthreads)

    from mrs_pipe.workflows.workflows import get_hermes_proc_wf

    basic_proc = get_hermes_proc_wf(args.subject, args.bids_dir)
    basic_proc.run()