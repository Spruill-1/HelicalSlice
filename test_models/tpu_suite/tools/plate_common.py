"""Where a plate builder is allowed to write.

A builder gets run for two different reasons: to VERIFY a code change - slice it, measure it,
read the numbers - and to PUBLISH a file that is going to be printed. Both used to write
straight to sliced/, so verifying a change silently overwrote the artifact that matched the
plastic on the bench. That happened three times in this campaign. The third cost a detour
through git to recover the exact file that had actually been printed, and in between the
repo's ledger quietly disagreed with the parts in hand.

Verifying is now the default and writes to a scratch directory git ignores. Touching sliced/
requires `--publish`, which is a deliberate act meaning "this is the file I intend to print".
"""
import os
import sys


def resolve_outdir(suite, argv=None):
    """Return (output_dir, publishing). Consumes `--publish` from argv so the caller's own
    argument parsing does not have to know about it."""
    argv = sys.argv if argv is None else argv
    if "--publish" in argv:
        argv.remove("--publish")
        out = os.path.join(suite, "sliced")
        print(f"[publish] writing the printable artifact to {out}")
        return out, True
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_verify")
    os.makedirs(out, exist_ok=True)
    print(f"[verify]  writing to {out} - sliced/ untouched (pass --publish to update it)")
    return out, False


def inject_env(outdir):
    """Environment for inject_thumbnails.py so it finds the file wherever it was written."""
    return dict(os.environ, TPU_SLICED_DIR=os.path.abspath(outdir))
