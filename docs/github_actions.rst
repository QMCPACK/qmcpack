.. _github_actions:

Github Actions, QMCPACK, and You!
=================================

This guide will cover the purpose and usual interactions a QMCPACK contributor might experience in regards to the Github Actions.  For how Github Actions accomplishes these things and other implementation details please refer to the offical `Github Actions Docs <https://docs.github.com/en/actions/guides>`_ and our scripts, which at the time of writing this doc are located `here <https://github.com/QMCPACK/qmcpack/tree/develop/tests/test_automation/github-actions/ci>`_.

********************************
So, What Even is Github Actions?
********************************

Good question! Github Actions is an event driven automation tool that allows us to automatically execute commands in response to QMCPACK repo related actions.  For example, merging a branch into master might then trigger our test scripts to run.

**********************************
Neat! How We are Using It and Why?
**********************************

One of the biggest hows and whys is that it saves time! There are certain jobs we always trigger upon merge, and having a designated person trigger and manage them takes a lot of time off their hands.  Time they could be spending actually improving the project with. 

Currently we are using Github Actions to automatically handle the following jobs:

* TODO: Add a list of jobs here
* Configure qmcpack using cmake out-of-source builds 
* Sanitize with clang compilers
* Build using ninja 
* Run deterministic tests
* Generate coverage reports
* Coming Soon: Install the library

detail each job

*****************************************
This is Pretty Cool, How Do I Contribute?
*****************************************

TODO: Firgure out examples where contributors might wanna add their own jobs and stuff, and how exactly they're supposed to do that.
TODO: Maybe layout some standards to keep everything clean and managable?
TODO: Review process for contributions? (security and such?)
TODO: Are we going to cover the different external runners in this and how to access them in the CI?
