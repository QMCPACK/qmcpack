#! /usr/bin/env -S uv run --script

import sys
from datetime import datetime
from pathlib import Path


def write_machine_table(output_file: Path):
    from nexus.machines import Machine, Supercomputer

    all_scs: list[Supercomputer] = []
    for val in Machine.machines.values():
        if not isinstance(val, Supercomputer):
            continue

        all_scs.append(val)

    output_text = (
        f".. Generated with the Python script at 'nexus/dev_utils/{Path(__file__).name}'\n\n"
        ".. currentmodule:: nexus.machines\n\n"
        ".. note::\n\n"
        f"    Table last updated on {datetime.now().strftime('%B %d, %Y')}\n\n"
    )

    output_text += (
        "+"+(28*"-")+
        "+"+( 7*"-")+
        "+"+( 9*"-")+
        "+"+(18*"-")+
        "+"+(14*"-")+
        "+"+(12*"-")+
        "+"+(14*"-")+
        "+"+(15*"-")+
        "+"+(13*"-")+"+\n"
    )
    output_text += (
        f"| {'Name':^26} "
        f"| {'Nodes':^5} "
        f"| {'Sockets':^7} "
        f"| {'Cores per Socket':^16} "
        f"| {'RAM per Node':^12} "
        f"| {'Queue Size':^10} "
        f"| {'Job Launcher':^12} "
        f"| {'Status Queuer':^13} "
        f"| {'Job Remover':^11} |\n"
    )
    output_text += (
        "+"+(28*"=")+
        "+"+( 7*"=")+
        "+"+( 9*"=")+
        "+"+(18*"=")+
        "+"+(14*"=")+
        "+"+(12*"=")+
        "+"+(14*"=")+
        "+"+(15*"=")+
        "+"+(13*"=")+"+\n"
    )

    for sc in all_scs:
        output_text += (
            f"| {f':py:class:`~.{type(sc).__name__}`':<26} "
            f"| {sc.nodes:>5} "
            f"| {sc.sockets_per_node:^7} "
            f"|      {sc.cores_per_socket:>4}        "
            f"|  {sc.ram_per_node:>6}      "
            f"|  {sc.queue_size:>6}    "
            f"| {sc.sub_launcher:^12} "
            f"| {sc.queue_querier:^13} "
            f"| {sc.job_remover:^11} |\n"
        )
        output_text += (
            "+"+(28*"-")+
            "+"+( 7*"-")+
            "+"+( 9*"-")+
            "+"+(18*"-")+
            "+"+(14*"-")+
            "+"+(12*"-")+
            "+"+(14*"-")+
            "+"+(15*"-")+
            "+"+(13*"-")+"+\n"
        )

    output_file.write_text(output_text)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        print("Useage: `./write_machine_table.py`")
        sys.exit(1)

    output_file = Path(__file__+"/../../docs/_static/machine_table.rst").resolve()
    write_machine_table(output_file)
    print(f'Updated machines table written to {output_file}')
