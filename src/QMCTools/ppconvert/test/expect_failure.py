#!/usr/bin/env python3

"""Run a command that must fail normally and emit an expected message.

This is stricter than CTest's WILL_FAIL property because it also rejects signal
termination and requires a specific diagnostic in the command output.

Example:
    python3 expect_failure.py --contains "Invalid option" -- ./program --bad-option
"""

import argparse
import subprocess
import sys


def main():
    """Run a command and require a normal failure containing expected output."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--contains", required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()

    command = (
        args.command[1:] if args.command and args.command[0] == "--" else args.command
    )
    if not command:
        parser.error("a command is required after --")

    result = subprocess.run(command, capture_output=True, text=True)
    output = result.stdout + result.stderr
    sys.stdout.write(output)

    if result.returncode == 0:
        print("Command unexpectedly succeeded", file=sys.stderr)
        return 1
    if result.returncode < 0:
        print(f"Command terminated by signal {-result.returncode}", file=sys.stderr)
        return 1
    if args.contains not in output:
        print(f"Expected output not found: {args.contains}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
