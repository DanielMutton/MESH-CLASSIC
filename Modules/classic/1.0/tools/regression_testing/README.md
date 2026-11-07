Regression Testing Toolkit {#regTest}
========

This collection of scripts works together to perform regression testing.

# What is Regression Testing?

Regression testing is a form of software testing that entails re-running tests on software after changes are introduced. The goal is to ensure that the same tests that passed before the changes, still pass after the changes. There are many different types of regression tests, but this particular suite is checksum-based. This means that any change to the logic of CLASSIC will radically change the outcome of the tests.

# What do these tests tell us?

The `verification.sh` script compares the checksum produced by CLASSIC to the checksum from a previous run. If they are identical, it means the logic of the program is unchanged. This is extremely valuable information when we only want to change the structure of the program, such as refactoring, commenting, and adding modules that are not yet being called.

If a developer is adding to or changing the logic of the program, these regression tests will fail, and that is to be expected. Once the changes are fully implemented, new checksums can be generated to compare future runs to.

# How to run the tests

This testing suite was designed to be run automatically in the continuous integration pipeline implemented on the supercomputer. Manual checksum testing can be done, but will require significant effort from the tester.

## Running the tests on the supercomputer

The regression tests will be performed whenever a commit is pushed to the code repository. If you are using your own fork of the repository, it is crucial that pipelines are enabled and that a runner is available to execute the commands. As of 2019/07/01, Shared Services has an account with a shared runner to be used for in-house projects, including CLASSIC. This is the runner that is enabled on the main CLASSIC repository. It is recommended that this same runner is used for other forks.

When you are finished with your changes to the code, commit and push them. The pipeline status can then be viewed from the GitLab web interface. If the pipeline passes, the regression tests were successful (and thus the checksums were unchanged). Again, if the changes to the code were meant to change the internal logic, then the checksums will fail, and that is ok. However, if the goal of the changes was structural, and the testing suite fails, then the changes should be reexamined.

## Running the tests manually

If these tests are being run outside of ECCC's supercomputer environment, then the only necessary file is `regtest.py`. A baseline checksum needs to be established, so this is best done before any changes are introduced. Run CLASSIC like normal, and note the output directory. Open `regtest.py`, and change the `output_directory` variable to the absolute path of your output directory. Then run `python3 regtest.py`. If you navigate to the output directory, there should now be a `checksums.csv` file. Copy this file somewhere else for reference.

After you have made non-logical changes to the code, run CLASSIC again. Be sure to use the same job options file, with the same input files. Repeat the process of running `regtest.py` on the output. You should have a new `checksums.csv` file in your output directory. The only step left is to use `diff [checksum_file_1] [checksum_file_2]` with the before-and-after checksum files. An empty diff result indicates a pass of the regression test (identical checksum files).

# Creating new test cases (supercomputer)

At the time of writing, there are currently three test cases being run in the pipeline: ctem_on, ctem_off, and ctem_compete. Each of these references a particular job options file as well as tasks in the pipeline. In order to create new test cases, we will add to these.

* First, a job options file should be created. The current runpath for checksum testing is `/space/hall1/sitestore/eccc/crd/ccrp/scrd530/classic_checksums/`. Navigate to this runpath, and into the `input_files` subdirectory. For simplicity's sake, it is best to copy one of the existing job options files (to keep the same input netcdf files) and simply tune the parameters to the new test case. The new file should be named `job_options_[name_of_test_case].txt`.
* Next, we will make a new submission script. Navigate to `tools/regression_testing/submission_scripts/` in your repository, and copy one of the existing submission scripts. Rename it to `job_submit_[name_of_test_case].sh`. Only two lines in the submission script must be modified, and both are at the top:
 * `/space/hall1/sitestore/eccc/crd/ccrp/scrd530/classic_checksums/[name_of_test_case]`
 * `/space/hall1/sitestore/eccc/crd/ccrp/scrd530/classic_checksums/input_files/job_options_[name_of_test_case].txt`
* Finally, we will create a new task by modifying `.gitlab-ci.yml` in the root of the repository. Again, this is best done by copying and modifying an existing test case. Two tasks in total need to be added: a run and a test. As an example, copy the `run:ctem_on` and `regression_test:ctem_on` jobs. Wherever you see `ctem_on` in these new jobs, change it to `[name_of_test_case]`. Push the changes and see if your new test case works!
