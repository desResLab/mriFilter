#include "mriProgramOptions.h"
#include "mriConstants.h"

#include <stdlib.h>
#include <getopt.h>

mriProgramOptions::mriProgramOptions(){
}


int mriProgramOptions::getCommadLineOptions(int argc, char **argv){
  int c;
  while (1){
    static struct option long_options[] =
    {
      {"pressure",  no_argument,       0, 1},
      {"template",  no_argument,       0, 2},
      {"input",     required_argument, 0, 'i'},
      {"output",    required_argument, 0, 'o'},
      {"smpfilter", required_argument, 0, 3},
      {"stats",     required_argument, 0, 4},
      {"matrix",    required_argument, 0, 5},
      {"random",    required_argument, 0, 6},
      {"crop",      required_argument, 0, 7},
      {"stream1",   required_argument, 0, 8},
      {"stream2",   required_argument, 0, 9},
      {"threshold", required_argument, 0, 10},
      {"reynolds",  required_argument, 0, 11},
      {"coeffs",    required_argument, 0, 12},
      {"gradient",  required_argument, 0, 13},
      {"vortex",    required_argument, 0, 14},
      {"writexp",   required_argument, 0, 15},
      {"test",      no_argument,       0, 16},
      {0, 0, 0, 0}
    };

    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "abc:d:f:",
                     long_options, &option_index);

    /* Detect the end of the options.*/
    if (c == -1){
      break;
    }

    // Check which option was selected
    switch (c){
      case 'i':
        // Read the input file name

        break;

      case 0:
        /* If this option set a flag, do nothing else now. */
        if (long_options[option_index].flag != 0){
          break;
        }
        printf ("option %s", long_options[option_index].name);
        if (optarg){
          printf (" with arg %s", optarg);
        }
        printf ("\n");
        break;
      case 'a':
        puts ("option -a\n");
        break;
      case 'b':
        puts ("option -b\n");
        break;
      case 'c':
        printf ("option -c with value `%s'\n", optarg);
        break;
      case 'd':
        printf ("option -d with value `%s'\n", optarg);
        break;
      case 'f':
        printf ("option -f with value `%s'\n", optarg);
        break;
      case '?':
        /* getopt_long already printed an error message. */
        break;
      case 16:
        runMode = rmSOLVEPOISSON;
        break;
      default:
        abort ();
      }
    }
    /* Print any remaining command line arguments (not options). */
    /*if (optind < argc){
      printf ("non-option ARGV-elements: ");
      while (optind < argc){
        printf ("%s ", argv[optind++]);
      }
      putchar ('\n');
    }*/
  return 0;
}
