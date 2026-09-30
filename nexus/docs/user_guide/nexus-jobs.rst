.. _nexus-jobs:

Detailed Nexus Job Usage
========================

This note describes how to use the ``job`` function in Nexus scripts. User scripts import `job` as follows:

.. code-block:: python

   from nexus import settings, job, run_project


Setting the Default Machine
---------------------------

Most example scripts set a global machine once with ``settings``, then let each
job inherit it:

.. code-block:: python

   settings(
       results = '',
       machine = 'ws16',
       )


Common workstation choices are ``ws1``, ``ws4``, ``ws8``, and ``ws16``. You may also use``machine='ws'`` or
``machine='workstation'`` to auto-detect physical CPU cores and register a matching workstation.

More commonly, a named supercomputer is used a specified in the machine entry.
For cluster or supercomputing machines that require an allocation, set the account globally:

.. code-block:: python

   settings(
       machine = 'frontier',
       account = 'ABC123',
       )


You can override the global machine or account on a single job:

.. code-block:: python

   job(machine='polaris', account='ABC123', nodes=2, app='qmcpack')


Attaching Jobs to Simulations
-----------------------------

Jobs are normally passed to ``generate_*`` functions through the ``job`` keyword:

.. code-block:: python

   scf = generate_pwscf(
       identifier = 'scf',
       path       = 'diamond/scf',
       job        = job(cores=16, app='pw.x'),
       )


When the simulation object is created, Nexus copies the job, assigns output
files and directories, fills in missing application commands, and asks the
machine to finalize core/node/process details.

Because the job is copied, reusable job templates are fine and appear often in
the examples:

.. code-block:: python

   scf_job = job(cores=16, threads=16)
   qmc_job = job(cores=16, threads=4, app='qmcpack')
   
   scf = generate_quantum_package(
       identifier = 'scf',
       path       = 'h2o/scf',
       job        = scf_job,
       )
   
   vmc = generate_qmcpack(
       identifier   = 'vmc',
       path         = 'h2o/vmc',
       job          = qmc_job,
       dependencies = (scf, 'orbitals'),
       )


Choose Cores, Nodes, and Threads
--------------------------------

For workstation examples, use ``cores``:

.. code-block:: python

   job(cores=16)


On a workstation, Nexus turns this into MPI ranks through the machine launcher.
For the default workstation launcher, ``job(cores=16)`` runs like
``mpirun -np 16 ...``.

For batch machines, use ``nodes`` when you want whole nodes:

.. code-block:: python

   job(nodes=2, hours=6, minutes=30)


``Supercomputer.process_job`` fills in the remaining resource fields from the
machine definition. If you specify ``nodes``, it computes total cores from the
machine's cores per node. If you specify ``cores``, it computes the node count
needed to cover those cores.

Use ``threads`` for OpenMP or other threaded codes:

.. code-block:: python

   job(cores=16, threads=4, app='qmcpack')


This requests 16 total cores with 4 threads per process. Nexus sets
``OMP_NUM_THREADS`` and computes the number of processes as ``cores / threads``
where possible. The examples use this pattern for QMCPACK hybrid jobs.

MPI Process Layout
^^^^^^^^^^^^^^^^^^

Nexus uses three related fields to describe the finalized MPI layout:

- ``processes`` is the total number of MPI ranks in the job. On a
  supercomputer, Nexus normally derives it as ``cores / threads``.
- ``processes_per_node`` is the MPI ranks placed on each allocated node. This
  is the main field to set when the default per-node layout is not appropriate.
  If it is set, Nexus uses ``nodes * processes_per_node`` as the total process
  count.
- ``processes_per_proc`` is the MPI ranks per physical processor package,
  usually a socket. Nexus derives it only when the ranks divide evenly across
  all processor packages; otherwise it is ``None``. It is mainly a description
  of the final layout, although some ``aprun`` machines use it to form ``-S``.

For example, use fewer MPI ranks per node with more OpenMP threads by setting
``processes_per_node`` explicitly:

.. code-block:: python

   job(nodes=2, threads=4, processes_per_node=8, app='qmcpack')


This creates 16 MPI ranks total, with four threads per rank. Make sure that
``processes_per_node * threads`` fits the CPU cores available on a node and
obeys any machine-specific CPU or GPU limits. In most cases, specify ``cores``
or ``nodes`` with ``threads`` and let Nexus derive the process layout; use
``processes_per_node`` only to select a deliberate rank/thread decomposition.

When using complex node, thread and process settings, it is highly recommended to verify that the generated job script is correct.

Running Serial or Threaded-Serial (non-MPI) Codes
-------------------------------------------------

Use ``serial=True`` when an application should not be launched with MPI:

.. code-block:: python

   job(serial=True, app='python3')


This is common for generic Python/Bash workflows and PySCF examples. With
``serial=True``, Nexus suppresses the MPI launcher when there is one process, so
the command is just the application command.

Threaded serial jobs are also supported:

.. code-block:: python

   job(serial=True, threads=16)


The PySCF examples use this form because PySCF do not support MPI by default, but is
highly threaded. In this case Nexus still exports ``OMP_NUM_THREADS=16``.

Specifying the Application
--------------------------

Use ``app`` for the executable name or path:

.. code-block:: python

   job(cores=16, app='pw.x')
   job(cores=1,  app='pw2qmcpack.x')
   job(cores=16, app='qmcpack')
   job(cores=8,  app='rmg-cpu')


Many simulation generators can supply a default executable command. In that
case the examples sometimes omit ``app``:

.. code-block:: python

   job(cores=16)


Use explicit ``app`` when you need a non-default executable, a converter, or a
full path:

.. code-block:: python

   job(cores=16, app='/path/to/qe/bin/pw.x')


Internally, ``app`` sets ``app_name``. During simulation initialization, Nexus
combines the simulation input file with the application name to form
``app_command``. Supplying ``app_command`` directly via `job` overrides the command generated by Nexus, substituting in the supplied command.  Use with caution.

Setting Job Walltimes
---------------------

Use ``hours``, ``minutes``, ``seconds``, and ``days``:

.. code-block:: python

   job(nodes=1, app='pw.x', hours=1)
   job(nodes=2, hours=6, minutes=30)


Nexus normalizes time values. For example, ``minutes=127`` becomes 2 hours and 7
minutes. Machine headers convert this time into the scheduler-specific format
such as PBS, Slurm, LoadLeveler, or LSF walltime.

Adding Environment Variables and Shell Snippets to Jobs
-------------------------------------------------------

Use ``env`` for environment variables that should be available to the job:

.. code-block:: python

   job(
       cores = 16,
       app   = 'qmcpack',
       env   = dict(QMC_DATA='data'),
       )


Nexus automatically sets ``OMP_NUM_THREADS`` from ``threads``.
For supercomputers, environment variables are written as exports in the batch
script. For workstations, they are passed to the subprocess environment.

Use ``presub`` and ``postsub`` for shell text immediately before or after the run
command:

.. code-block:: python

   job(
       cores  = 16,
       app    = 'qmcpack',
       presub = 'module load qmcpack',
       )


For supercomputers, these lines are placed in the batch script. For
workstations, they are included in the command block run locally.

Common uses for these options are loading modules or activating a special Python environment for a particular job.

Passing Options to Jobs
-----------------------

Nexus tracks three option groups:

- ``app_options``: appended after the executable command, such as `pw.x` or `qmcpack`.
- ``run_options``: appended after the run launcher, such as
  ``mpirun``, ``srun``, or ``aprun``.
- ``sub_options``: appended after the submit launcher, such as ``qsub`` or
  ``sbatch``.

Examples usually let the machine set launcher options. Reach for these only
when you need a code-specific flag or a machine-specific override:

.. code-block:: python

   job(
       cores       = 16,
       app         = 'qmcpack',
       app_options = '--verbosity=high',
       )


Option inputs can be strings, lists of strings, or dictionaries. Dictionary
keys are used only for stable storage and ordering; the values should contain
the command-line text to write.

Using Batch-Specific Keywords When Needed
-----------------------------------------

Machine classes interpret batch fields in their own ``write_job_header`` or
processing hooks. Common fields include:

- ``queue``: queue or partition name.
- ``account``: allocation or project account.
- ``qos``: scheduler quality of service.
- ``constraint``: architecture or node-type constraint.
- ``gpus``: GPUs requested per node or per machine-specific convention.
- ``filesystems``: PBS filesystem request on machines that require it.
- ``reservation``: reservation name for machines that require reservations.

For example:

.. code-block:: python

   job(
       machine = 'frontier',
       account = 'ABC123',
       nodes   = 2,
       queue   = 'batch',
       hours   = 1,
       app     = 'qmcpack',
       )


If a machine has defaults for ``queue`` or similar fields, leave them unset
unless you need something specific.

Local Execution and Skip-Submit
-------------------------------

``local=True`` tells a simulation to execute immediately instead of entering the
machine queue:

.. code-block:: python

   job(cores=4, local=True)


This is different from ``serial=True``. ``serial=True`` controls whether an MPI
launcher is used. ``local=True`` controls whether Nexus queues/submits the job or
executes it directly from the simulation.

The examples more often use project-level controls such as ``generate_only``,
``status_only``, and ``skip_submit`` in ``settings`` or command-line flags. These are
workflow controls, not job-resource settings.

Using Full Commands to Override Job Defaults
--------------------------------------------

``full_command`` bypasses the normal launcher/application construction:

.. code-block:: python

   job(full_command='custom_launcher --flag input.in > out 2> err')


Only use this when Nexus cannot express the command with ``app``, ``app_options``, and ``run_options``. When ``full_command`` is
set, Nexus does not add the machine launcher or application options for you. Use of the ``full_command`` option is discouraged
because it reduces the portability of scripts.

Common Recipes from the Examples
--------------------------------

Workstation MPI code (16 MPI tasks only):

.. code-block:: python

   job(cores=16, app='pw.x')


Converter or small helper:

.. code-block:: python

   job(cores=1, app='pw2qmcpack.x')


Hybrid MPI/OpenMP QMCPACK (4 MPI tasks each with 4 OpenMP threads):

.. code-block:: python

   job(cores=16, threads=4, app='qmcpack')


Threaded PySCF without MPI:

.. code-block:: python

   job(serial=True, threads=16)


Generic Python or Bash task:

.. code-block:: python

   job(serial=True, app='python3')
   job(serial=True, app='bash')


Reusable templates:

.. code-block:: python

   qp_job  = job(cores=16, threads=16)
   c4q_job = job(cores=1)
   qmc_job = job(cores=16, threads=16)


Batch job with scheduler metadata:

.. code-block:: python

   job(
       machine = 'polaris',
       account = 'ABC123',
       nodes   = 2,
       queue   = 'prod',
       hours   = 2,
       app     = 'qmcpack',
       )


Verifying and Debugging Generated Scripts
-----------------------------------------

When creating new Jobs or using new machine, it is highly recommended to verify that the created job scripts are correct. As with
any job, incorrect settings can lead to very slow or failed runs and wasted resources.

You can run a script with generation mode and inspect the generated submission files manually:

.. code-block:: bash

   ./my_nexus_script.py --generate_only


Alternatively, within Python, you can print out the submission file contents ahead of time:

.. code-block:: python

 scf = generate_pwscf(
       identifier = 'scf',
       path       = 'diamond/scf',
       job        = job(machine = 'polaris', 
                        nodes   = 2, 
                        hours   = 1, 
                        app     = 'pw.x'),
       )

  print(scf.job.write())
  exit()

