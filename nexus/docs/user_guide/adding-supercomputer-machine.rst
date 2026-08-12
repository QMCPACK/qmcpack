.. _adding-supercomputer-machine:

Adding a Supercomputer Machine
==============================

This section describes how to add a permanent machine definition to ``machines.py``. These definitions are used by Nexus to generate
job submission scripts, and subsequently to help submit and monitor jobs. Within Nexus, ``Supercomputer`` refers to any machine that
uses a batch scheduler such at Slurm or PBS and not the size of the machine. This could include even a small workstation where a
batch scheduler is installed and used to schedule jobs overnight. If a new machine only requires common or default behaviors, the
new machine definition will require minimal Python code due to the functionality provided by the ``Supercomputer`` base class. Runs
on machines that do not have a batch scheduler should use one of the workstation machine definitions.

We encourage submission of new permanent machine definitions to the Nexus repository on GitHub to reduce the need for ongoing local
maintenance and to share them with other users. For a one-off local computational cluster, the same class and instantiation patterns
can be used in``~/.nexus/local_machines.py``. 

In the following we describe a step-by-step process for adding a new machine. 

Choose the Closest Existing Machine
-----------------------------------

To minimize the work required, start by finding a machine that uses the same scheduler and run launcher.
For example:

- Slurm plus ``srun``: see ``Andes``, ``Rhea``, ``Frontier``, ``Leonardo``.
- Slurm plus ``mpirun``: see ``CadesSlurm``, ``Tomcat3``, ``Improv``.
- PBS plus ``mpiexec``: see ``Polaris``, ``Aurora``.
- PBS plus ``aprun``: see ``BlueWatersXE``, ``BlueWatersXK``, ``Theta``.
- LSF plus ``jsrun`` or ``lrun``: see ``Summit``, ``Lassen``.

If the header and launcher behavior are shared by several machines, consider
deriving from the existing intermediate class, such as ``NerscMachine`` or
``SnlMachine``. Otherwise, derive directly from ``Supercomputer``.

Adding the Machine Class
------------------------

Place the new class near similar machines in ``machines.py``. At minimum, define
a unique lower-case ``name``, the account/capability flags, and
``write_job_header``.

.. code-block:: python

   class NewMachine(Supercomputer):
       name = 'newmachine'
       requires_account = True
       batch_capable    = True
   
       def post_process_job(self,job):
           if len(job.run_options)==0:
               job.run_options.add(
                   N = '-N {}'.format(job.nodes),
                   n = '-n {}'.format(job.processes),
                   c = '-c {}'.format(job.threads),
                   )
           #end if
       #end def post_process_job
   
       def write_job_header(self,job):
           if job.queue is None:
               job.queue = 'regular'
           #end if
           c  = '#!/bin/bash\n'
           c += '#SBATCH -A {}\n'.format(job.account)
           c += '#SBATCH -p {}\n'.format(job.queue)
           c += '#SBATCH -J {}\n'.format(job.name)
           c += '#SBATCH -t {}\n'.format(job.sbatch_walltime())
           c += '#SBATCH -N {}\n'.format(job.nodes)
           c += '#SBATCH --ntasks-per-node={}\n'.format(job.processes_per_node)
           c += '#SBATCH --cpus-per-task={}\n'.format(job.threads)
           c += '#SBATCH -o {}\n'.format(job.outfile)
           c += '#SBATCH -e {}\n'.format(job.errfile)
           if job.user_env:
               c += '#SBATCH --export=ALL\n'
           else:
               c += '#SBATCH --export=NONE\n'
           #end if
           return c
       #end def write_job_header
   #end class NewMachine

Please note: the only required class method is ``write_job_header``.  Defining the ``post_process_job`` function can be useful in special cases, see below.

Important class details:

- ``name`` is the key used by ``Machine.get``, ``job(machine=...)``, and the tests.
  It must be unique in ``Machine.machines``.
- ``requires_account = True`` makes jobs require ``account``; the machine tests
  supply ``ABC123`` automatically for such machines.
- ``write_job_header`` returns only the batch-script header and any setup lines.
  ``Supercomputer.write_job`` appends environment exports and the run command.
- Use ``pre_process_job`` to set defaults or hardware variants before generic
  node/core calculations. Use ``post_process_job`` to add launcher options after
  ``job.nodes``, ``job.processes``, ``job.processes_per_node``, and ``job.threads``
  have been finalized.
- Keep ``process_job`` behavior idempotent. The tests call it more than once on
  already-processed jobs, so avoid appending duplicate options or mutating
  machine-wide state in a way that changes later jobs unexpectedly.

Understanding the role of the Base Class
----------------------------------------

``Supercomputer.process_job`` fills in missing ``cores`` or ``nodes``, computes
``processes``, ``processes_per_node``, ``processes_per_proc``, ``ppn``, applies the
machine account default, sets ``OMP_NUM_THREADS``, and then calls
``process_job_options``.

Default launcher handling is limited:

- ``mpirun`` adds ``-np <processes>``.
- ``mpiexec`` adds ``-n <processes>``.
- ``aprun`` adds ``-n <processes>`` and, for threaded jobs, ``-d <threads>``.
- ``runjob`` adds Blue Gene style ``--np``, ``-p``, ``$LOCARGS``, and ``--envs``.
- ``ibrun`` adds ``-n <processes> -o 0``.
- ``srun``, ``jsrun``, and ``lrun`` intentionally add nothing by default.

For launchers that need machine-specific options, override ``post_process_job``
or ``process_job_options``. Prefer ``post_process_job`` when you only need to add
or adjust ``job.run_options``; override ``process_job_options`` when the base
launcher behavior is not appropriate at all.

``Options.write()`` sorts option keys before building the command, so choose
stable keys if test output ordering matters.

Registering the Machine
-----------------------

At the bottom of ``machines.py``, add an instance to the ``#Known machines`` block:

.. code-block:: python

   #            nodes sockets cores ram qslots qlaunch qsubmit qstatus qdelete
   NewMachine(  1000,   2,    32, 256, 1000,  'srun', 'sbatch','squeue','scancel')


The positional arguments are:

- ``nodes``: total compute nodes.
- ``procs_per_node``: usually sockets or processor packages per node.
- ``cores_per_proc``: cores per socket/package.
- ``ram_per_node``: memory per node, in GB.
- ``queue_size``: maximum queued/running jobs Nexus should track.
- ``app_launcher``: command used in ``Job.run_command``.
- ``sub_launcher``: command used by ``sub_command`` and submission file names.
- ``queue_querier``: queue-status command/parser. Supported values include
  ``qstat``, ``qstata``, ``squeue``, ``sacct``, ``llq``, ``bjobs``, and ``test_query``.
- ``job_remover``: cancellation command. Existing support includes ``qdel`` and
  ``scancel``; add support if the scheduler needs another remover.

``cores_per_node`` is computed as ``procs_per_node * cores_per_proc``. Use the
actual socket/package layout when it matters for binding or
``processes_per_proc``.

Updating the Machine Tests
--------------------------

Adding the instance registers the machine in ``Machine.machines``, so the generic
machine tests will include it automatically.

Update ``tests/test_machines.py::test_job_run_command`` by adding entries to
``job_run_ref`` for the new machine. Unless the test has a special case for the
machine, add all six standard job shapes:

- ``n1``
- ``n1_p1``
- ``n2``
- ``n2_t2``
- ``n2_t2_e``
- ``n2_t2_p2``

For a permanent machine, also update ``test_write_job`` by adding the expected
batch script to ``job_write_ref``; otherwise that test will hit a missing
reference for the newly registered machine.

To print fresh references, run this command from the Nexus repository root:

.. code-block:: bash

   dev_utils/write-job-ref


This writes ``updated_job_ref_table.txt`` to the current working directory. The
file separates the entries for ``job_run_ref`` in ``test_job_run_command`` from the
entries for ``job_write_ref`` in ``test_write_job``. Copy only the lines for the new
machine into the corresponding table, then review them manually. The generated
commands and submission files are useful starting points, not a substitute for
checking scheduler syntax against the machine documentation.

Running the Tests
-----------------

From the repository root, run the machine tests that cover the new definition:

.. code-block:: bash

   python -m pytest \
     nexus/tests/test_machines.py::test_process_job \
     nexus/tests/test_machines.py::test_job_run_command \
     nexus/tests/test_machines.py::test_write_job

If the new machine requires an extra mandatory job field, either provide a safe
default in the machine class or add a narrow special case in the relevant test,
as ``summit`` and ``flight`` already do.

If you are also using and configuring QMCPACK alongside Nexus, you can more simply run the machines test using ``ctest -R ntest_nexus_machines``, or all the Nexus tests with ``ctest -R nexus``,
from the QMCPACK build directory.
