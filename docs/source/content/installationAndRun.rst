Installing and Running smpFilter
================================

Installation
^^^^^^^^^^^^
First, clone the git repository in a local folder:

  git clone git@bitbucket.org:danschi/mpfilter.git 

Create a directory for out-of-source build (e.g., smpBin):

  mkdir smpBin

Enter this folder 

  cd smpBin

run cmake with the following command:

  cmake ../mpfilter/

Run make:

  make (or make -jn for parallel compilation with n processes)

#### Running the code

The preferred way to run the code is by reading a command file:

  mpFilterApp -c command_file_name

