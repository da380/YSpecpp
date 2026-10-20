"""The command `yspec-plot`."""

import argparse
import sys

from .runs import SeismogramRun, SnapshotRun

USAGE = """Pictures of what the YSpec programs write, from the input file of a run.

  yspec-plot section INPUT [--component Z] [--common] [--gain G]
  yspec-plot traces INPUT [--receivers 1 2 ...] [--components ZRT]
  yspec-plot compare INPUT INPUT ... [--labels ...] [--receiver 1]
  yspec-plot snapshot INPUT [--component up] [--common]

section: a record section of one component; traces: the components at a
few receivers; compare: one receiver in several runs, with their
differences from the last; snapshot: a section through the earth, or a map
around the source. Components of seismograms are Z (up), N, E, R (radial,
away from the source) and T (transverse, 90 degrees anticlockwise from
radial seen from above). With several sources, --source chooses one by
name. Without -o the figure is shown.
"""


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="yspec-plot", description=USAGE,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)

    def command(name, help):
        sub = commands.add_parser(name, help=help)
        sub.add_argument("-o", "--output", help="the file to save it in")
        sub.add_argument("--source", help="the source, by name")
        return sub

    sub = command("section", "a record section")
    sub.add_argument("input")
    sub.add_argument("--component", default="Z")
    sub.add_argument("--common", action="store_true",
                     help="scale every trace alike")
    sub.add_argument("--gain", type=float, default=1,
                     help="amplify, clipping where traces would overlap")
    sub = command("traces", "the components at a few receivers")
    sub.add_argument("input")
    sub.add_argument("--receivers", type=int, nargs="+")
    sub.add_argument("--components", default="ZRT")
    sub = command("compare", "one receiver in several runs")
    sub.add_argument("inputs", nargs="+")
    sub.add_argument("--labels", nargs="+")
    sub.add_argument("--receiver", type=int, default=1)
    sub.add_argument("--component", default="Z")
    sub = command("snapshot", "a section through the earth or a map")
    sub.add_argument("input")
    sub.add_argument("--component", default="up",
                     choices=["up", "away", "transverse"])
    sub.add_argument("--common", action="store_true",
                     help="scale every time alike")
    sub.add_argument("--no-coastlines", action="store_true",
                     help="leave the coastlines off a map")
    args = parser.parse_args(argv)

    import matplotlib
    if args.output:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from .seismograms import compare, record_section, traces
    from .snapshots import snapshots
    try:
        if args.command == "section":
            fig = record_section(
                SeismogramRun(args.input, args.source), args.component,
                common=args.common, gain=args.gain)
        elif args.command == "traces":
            fig = traces(SeismogramRun(args.input, args.source),
                         args.receivers, args.components)
        elif args.command == "compare":
            fig = compare(
                [SeismogramRun(path, args.source) for path in args.inputs],
                args.labels, receiver=args.receiver,
                component=args.component)
        else:
            fig = snapshots(SnapshotRun(args.input, args.source),
                            args.component, common=args.common,
                            coastlines=not args.no_coastlines)
    except (OSError, ValueError, KeyError) as error:
        sys.exit(f"yspec-plot: {error}")
    if args.output:
        fig.savefig(args.output, dpi=150, bbox_inches="tight")
        print(f"wrote {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
