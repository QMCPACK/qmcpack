.. _adding-supercomputer-machine:

.. currentmodule:: nexus.machines

Adding a Supercomputer Machine
==============================

This section describes how to add a permanent machine definition to :py:mod:`~nexus.machines`. These definitions are used by Nexus to generate
job submission scripts, and subsequently to help submit and monitor jobs. Within Nexus, :py:class:`Supercomputer` refers to any machine that
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

Place the new class near similar machines in :py:mod:`~.nexus.machines`. At minimum, define
a unique lower-case :py:attr:`~Supercomputer.name`, the account/capability flags, and
``write_job_header`` method.

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

Please note: the only required class method is :py:meth:`~Supercomputer.write_job_header`.  Defining the :py:meth:`Supercomputer.post_process_job` function can be useful in special cases, see below.

Important class details:

- :py:attr:`~Supercomputer.name` is the key used by :py:meth:`~Machine.get`, ``job(machine=...)``, and the tests.
  It must be unique in :py:attr:`Machine.machines`.
- ``requires_account = True`` makes jobs require ``account``; the machine tests
  supply ``ABC123`` automatically for such machines.
- :py:meth:`~Supercomputer.write_job_header` returns only the batch-script header and any setup lines.
  :py:meth:`Supercomputer.write_job` appends environment exports and the run command.
- Use ``pre_process_job`` to set defaults or hardware variants before generic
  node/core calculations. Use :py:meth:`~Supercomputer.post_process_job` to add launcher options after
  ``job.nodes``, :py:attr:`Job.processes`, ``job.processes_per_node``, and ``job.threads``
  have been finalized.
- Keep :py:meth:`Supercomputer.process_job` behavior idempotent. The tests call it more than once on
  already-processed jobs, so avoid appending duplicate options or mutating
  machine-wide state in a way that changes later jobs unexpectedly.

Understanding the role of the Base Class
----------------------------------------

The :py:meth:`Supercomputer.process_job` method fills in missing :py:attr:`Job.cores` or :py:attr:`Job.nodes`, computes
:py:attr:`Job.processes`, :py:attr:`Job.processes_per_node`, :py:attr:`Job.processes_per_proc`, :py:attr:`Job.ppn`, applies the
machine account default, sets ``OMP_NUM_THREADS``, and then calls
:py:meth:`~Supercomputer.process_job_options`.

Default launcher handling is limited:

- ``mpirun`` adds ``-np <processes>``.
- ``mpiexec`` adds ``-n <processes>``.
- ``aprun`` adds ``-n <processes>`` and, for threaded jobs, ``-d <threads>``.
- ``runjob`` adds Blue Gene style ``--np``, ``-p``, ``$LOCARGS``, and ``--envs``.
- ``ibrun`` adds ``-n <processes> -o 0``.
- ``srun``, ``jsrun``, and ``lrun`` intentionally add nothing by default.

For launchers that need machine-specific options, override :py:meth:`Supercomputer.post_process_job`
or :py:meth:`Supercomputer.process_job_options`. Prefer :py:meth:`Supercomputer.post_process_job` when you only need to add
or adjust :py:attr:`Job.run_options`; override :py:meth:`Supercomputer.process_job_options` when the base
launcher behavior is not appropriate at all.

The :py:meth:`Options.write()` method sorts option keys before building the command, so choose
stable keys if test output ordering matters.

Registering the Machine
-----------------------

The :py:func:`register_supercomputer` decorator is used to automatically register a supercomputer in :py:attr:`Machine.machines`.

You must define the following parameters in the subclass you have created:

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
