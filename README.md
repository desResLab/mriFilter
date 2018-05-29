### smpFilter Application

This application implements the Sequential Matching Pursuit filter discussed in the following publication:

Schiavazzi D., Coletti F., Iaccarino G. and Eaton, J.K., **A matching pursuit approach to solenoidal filtering of three-dimensional velocity measurements**, _Journal of Computational Physics_, 263, 206--221, 2014.

### Installation prerequisites

#### Note for installation on the Certainty Cluster

The following modules should be loaded before running cmake and make to compile the application.

- **cmake/2.8.12.2** . The CMake utility.                 
- **mvapich2/2.1a-gcc-4.4.7** . MPI Libraries compiled with gcc.
- **intel/14** . The Intel compiler.
- **boost/1.55.0-mvapich2-2.0rc1-intel-14** . The Boost libraries.

### Installation

First, clone this git repository in a local folder:

Create a directory for out-of-source build (e.g., mriBin):

```
mkdir mriBin
```

Enter this folder 

```
cd mriBin
```

run cmake with the following command:

```
cmake ../mriFilter/
```

Run make:

```
make (or make -jn for parallel compilation with n processes)
```

#### Installing the documentation

The smpFilter documentation is written using Sphinx. To compile the documentation go to the docs folder:

```
cd docs
```

and run the make utility and specify the html target

```
make html
```

An **index.html** document will be created in the **docs/build/html** folder.

#### Running the code

The preferred way to run the code is by reading a command file:

```
mriFilter -c command_file_name
```
