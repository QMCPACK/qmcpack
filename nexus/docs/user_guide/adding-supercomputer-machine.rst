.. _adding-supercomputer-machine:

.. currentmodule:: nexus.machines

Adding a Supercomputer Machine
==============================

This section describes how to add a permanent machine definition to :py:mod:`~nexus.machines`. These definitions are used by Nexus to generate
job submission scripts, and subsequently to help submit and monitor jobs. Within Nexus, 'Supercomputer' refers to any machine that
uses a batch scheduler such at Slurm or PBS and not the size of the machine. This could include even a small workstation where a
batch scheduler is installed and used to schedule jobs overnight. If a new machine only requires common or default behaviors, the
new machine definition will require minimal Python code due to the functionality provided by the :py:class:`Supercomputer` base class. Runs
on machines that do not have a batch scheduler should use one of the workstation machine definitions.

We encourage submission of new permanent machine definitions to the Nexus repository on GitHub to reduce the need for ongoing local
maintenance and to share them with other users. For a one-off local computational cluster, the same class and instantiation patterns
can be used in ``~/.nexus/local_machines.py``.

In the following we describe a step-by-step process for adding a new machine. 


Choose the Closest Existing Machine
-----------------------------------

To minimize the work required, start by finding a machine that uses the same scheduler and run launcher.
For example:

- Slurm plus ``srun``: see :py:class:`Andes`, :py:class:`Rhea`, :py:class:`Frontier`, :py:class:`Leonardo`.
- Slurm plus ``mpirun``: see :py:class:`CadesSlurm`, :py:class:`Tomcat3`, :py:class:`Improv`.
- PBS plus ``mpiexec``: see :py:class:`Polaris`, :py:class:`Aurora`.
- PBS plus ``aprun``: see :py:class:`BlueWatersXE`, :py:class:`BlueWatersXK`, :py:class:`Theta`.
- LSF plus ``jsrun`` or ``lrun``: see :py:class:`Summit`, :py:class:`Lassen`.

If the header and launcher behavior are shared by several machines, consider
deriving from the existing intermediate class, such as :py:class:`NerscMachine` or
:py:class:`SnlMachine`. Otherwise, derive directly from :py:class:`Supercomputer`.


Adding the Machine Class
------------------------

Place the new class near similar machines in :py:mod:`~.nexus.machines`, and make sure
to add the :py:func:`register_supercomputer` decorator to the class.

Example Implementation
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    @register_supercomputer
    class NewMachine(Supercomputer):
        name = 'newmachine'
        requires_account = True
        batch_capable    = True
 
        nodes            = 3072
        sockets_per_node = 2
        cores_per_socket = 128
        ram_per_node     = 512
        queue_size       = 5000
        app_launcher     = "srun"
        sub_launcher     = "sbatch"
        queue_querier    = "squeue"
        job_remover      = "scancel"

        def post_process_job(self,job):
            if len(job.run_options)==0:
                job.run_options.add(
                    N = f'-N {job.nodes}',
                    n = f'-n {job.processes}',
                    c = f'-c {job.threads}',
                    )
            #end if
        #end def post_process_job

        def write_job_header(self,job):
            if job.queue is None:
                job.queue = 'regular'
            #end if
            c  = '#!/bin/bash\n'
            c += f'#SBATCH -A {job.account}\n'
            c += f'#SBATCH -p {job.queue}\n'
            c += f'#SBATCH -J {job.name}\n'
            c += f'#SBATCH -t {job.sbatch_walltime()}\n'
            c += f'#SBATCH -N {job.nodes}\n'
            c += f'#SBATCH --ntasks-per-node={job.processes_per_node}\n'
            c += f'#SBATCH --cpus-per-task={job.threads}\n'
            c += f'#SBATCH -o {job.outfile}\n'
            c += f'#SBATCH -e {job.errfile}\n'
            if job.user_env:
                c += '#SBATCH --export=ALL\n'
            else:
                c += '#SBATCH --export=NONE\n'
            #end if
            return c
        #end def write_job_header
    #end class NewMachine


Required Methods/Attributes
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automethod:: nexus.machines.Supercomputer.write_job_header
    :no-index:

.. autoattribute:: nexus.machines.Supercomputer.name
    :no-index:

.. autoattribute:: nexus.machines.Supercomputer.nodes
    :no-index:

.. autoattribute:: nexus.machines.Supercomputer.sockets_per_node
    :no-index:

.. autoattribute:: nexus.machines.Supercomputer.cores_per_socket
    :no-index:

.. autoattribute:: nexus.machines.Supercomputer.ram_per_node
    :no-index:

.. autoattribute:: nexus.machines.Supercomputer.queue_size
    :no-index:

.. autoattribute:: nexus.machines.Supercomputer.app_launcher
    :no-index:

.. autoattribute:: nexus.machines.Supercomputer.sub_launcher
    :no-index:

.. autoattribute:: nexus.machines.Supercomputer.queue_querier
    :no-index:

.. autoattribute:: nexus.machines.Supercomputer.job_remover
    :no-index:

.. autoattribute:: nexus.machines.Supercomputer.cores_per_node
    :no-index:


Optional Methods/Attributes
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automethod:: nexus.machines.Supercomputer.pre_process_job
    :no-index:

.. automethod:: nexus.machines.Supercomputer.post_process_job
    :no-index:

.. autoattribute:: nexus.machines.Machine.batch_capable
    :no-index:

.. autoattribute:: nexus.machines.Machine.requires_account
    :no-index:

.. autoattribute:: nexus.machines.Machine.executable_subfile
    :no-index:

.. autoattribute:: nexus.machines.Machine.redirect_output
    :no-index:

.. autoattribute:: nexus.machines.Machine.query_with_username
    :no-index:

.. autoattribute:: nexus.machines.Machine.special_bundling
    :no-index:

.. autoattribute:: nexus.machines.Machine.prefixed_output
    :no-index:

.. autoattribute:: nexus.machines.Machine.outfile_extension
    :no-index:

.. autoattribute:: nexus.machines.Machine.errfile_extension
    :no-index:

.. autoattribute:: nexus.machines.Machine.allow_warnings
    :no-index:

.. autoattribute:: nexus.machines.Machine.queue_configs
    :no-index:


.. note::

    The only required class method is :py:meth:`~Supercomputer.write_job_header`. Defining the :py:meth:`~Supercomputer.post_process_job` function can be useful in special cases, see below.

Important class details:

- :py:meth:`~Supercomputer.write_job_header` returns only the batch-script header and any setup lines.
  :py:meth:`~Supercomputer.write_job` appends environment exports and the run command.
- Keep :py:meth:`~Supercomputer.process_job` behavior idempotent. The tests call it more than once on
  already-processed jobs, so avoid appending duplicate options or mutating
  machine-wide state in a way that changes later jobs unexpectedly.


Understanding the role of the Base Class
----------------------------------------

The :py:meth:`~Supercomputer.process_job` method fills in missing :py:attr:`Job.cores` or :py:attr:`Job.nodes`, computes
:py:attr:`Job.processes`, :py:attr:`Job.processes_per_node`, :py:attr:`Job.processes_per_socket`, :py:attr:`Job.ppn`, applies the
machine account default, sets ``OMP_NUM_THREADS``, and then calls
:py:meth:`~Supercomputer.process_job_options`.

Default launcher handling is limited:

- ``mpirun`` adds ``-np <processes>``.
- ``mpiexec`` adds ``-n <processes>``.
- ``aprun`` adds ``-n <processes>`` and, for threaded jobs, ``-d <threads>``.
- ``runjob`` adds Blue Gene style ``--np``, ``-p``, ``$LOCARGS``, and ``--envs``.
- ``ibrun`` adds ``-n <processes> -o 0``.
- ``srun``, ``jsrun``, and ``lrun`` intentionally add nothing by default.

For launchers that need machine-specific options, override :py:meth:`~Supercomputer.post_process_job`
or :py:meth:`~Supercomputer.process_job_options`. Prefer :py:meth:`~Supercomputer.post_process_job` when you only need to add
or adjust :py:attr:`Job.run_options`; override :py:meth:`~Supercomputer.process_job_options` when the base
launcher behavior is not appropriate at all.

The :py:meth:`Options.write()` method sorts option keys before building the command, so choose
stable keys if test output ordering matters.


Updating the Machine Tests
--------------------------

Adding the :py:func:`register_supercomputer` decorator registers the machine in :py:attr:`Machine.machines`, so the generic machine tests will include it automatically.

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
as :py:class:`Summit` and :py:class:`Flight` already do.

If you are also using and configuring QMCPACK alongside Nexus, you can more simply run the machines test using ``ctest -R ntest_nexus_machines``, or all the Nexus tests with ``ctest -R nexus``,
from the QMCPACK build directory.
