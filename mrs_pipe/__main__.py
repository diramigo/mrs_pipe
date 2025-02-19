import argparse
import sys
from mrs_pipe.workflows import run_workflow

def main():
    parser = argparse.ArgumentParser(
        description="CLI para ejecutar scripts en mi_paquete."
    )
    parser.add_argument(
        "script", choices=["run_workflow"], help="El script a ejecutar"
    )
    parser.add_argument(
        "args", nargs=argparse.REMAINDER, help="Argumentos para el script seleccionado"
    )

    args = parser.parse_args()

    if args.script == "run_workflow":
        run_workflow.main(args.args)

if __name__ == "__main__":
    main()
